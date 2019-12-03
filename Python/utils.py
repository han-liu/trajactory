import os
import cv2
import glob
import pandas as pd
import numpy as np
import SimpleITK as sitk
import random
import matplotlib.pyplot as plt
from scipy.ndimage import zoom
from array2gif import write_gif
import cc3d


####################################
############## imgs ################
####################################

def readImg(imgFp):
    assert os.path.isfile(imgFp)
    return simg2img(sitk.ReadImage(imgFp))

def loadImg(imgFp, preprocess=False, targetSize=None):
    img = readImg(imgFp)
    if targetSize: img = resize(img, targetSize)
    if preprocess: img = preprocessImg(img)
    return img

def preprocessImg(img):
    img += 1024 # range -1024~600
    img = img.astype(dtype="uint16")
    return img

def saveImg(img, saveFp):
    sim = img2simg(img)
    sitk.WriteImage(sim, saveFp, True)
    print(f"Successfully saved image at {saveFp}")

def isotropic(imgFp, saveFp=None, zoomOrder=3):
    simg = sitk.ReadImage(imgFp)
    img = readImg(imgFp)
    spacing = simg.GetSpacing()
    shape = simg.GetSize()
    isoShape= tuple(map(lambda x, y: round(x*y+1e-5), spacing, shape))
    isoImg = resize(img, isoShape, zoomOrder)
    if saveFp: saveImg(isoImg, saveFp)   
    return isoImg

####################################
############## msks ################
####################################

def loadMsk(mskFp, preprocess=False, targetSize=None):
    msk = readImg(mskFp)
    if preprocess: msk = preprocessMsk(msk)
    msks = generateMsks(msk)
    if targetSize: msks = [resize(msk, targetSize, 0) for msk in msks]
    return msks

def preprocessMsk(msk):
    msk[np.where(msk!=1)]=0
    return msk

def generateMsks(msk):
    msks = []
    labels = np.unique(msk)
    for label in labels:
        x = msk.copy()
        if label != 0:
            x[np.where(x!=label)] = 0  
            x[np.where(x!=0)] = 1
        else:
            x[np.where(x!=0)] = 255
            x[np.where(x==0)] = 1
            x[np.where(x!=1)] = 0
        msks.append(x)
    return msks

####################################
############# utility ##############
####################################

def resize(img, targetSize, zoomOrder=3):
    return zoom(img, (targetSize[0]/img.shape[0],
                      targetSize[1]/img.shape[1],
                      targetSize[2]/img.shape[2]), order=zoomOrder)

def simg2img(simg):
    img = sitk.GetArrayFromImage(simg)
    img = np.swapaxes(img, 0, 2)
    img = np.swapaxes(img, 0, 1)
    return img

def img2simg(img):
    img = np.swapaxes(img, 0, 1)
    img = np.swapaxes(img, 0, 2)
    simg = sitk.GetImageFromArray(img)
    return simg

def get_id(imgFp):
    return os.path.splitext(os.path.basename(imgFp))[0]

def getLargestConnectedComp(msk):
    # msk: binary mask
    connectivity = 6  # only 26, 18, and 6 are allowed
    msk = cc3d.connected_components(msk.astype('uint8'), connectivity=connectivity)
    indices, counts = np.unique(msk, return_counts=True)
    indices = indices[1:]  # because zero is background
    counts = counts[1:]
    ccLabel = indices[np.argmax(counts)]
    msk[np.where(msk != ccLabel)] = 0
    msk[np.where(msk != 0)] = 1
    return msk

####################################
########### Create CSV #############
####################################

def generateCSV(imgDir, mskDir, saveDir):
    saveFp = os.path.join(saveDir, "img-msk.csv") 
    f = open(saveFp, "w")
    f.write("imgFp,mskFp\n")
    print("Generating csv file.....")
    imgFps = [os.path.join(imgDir, imgFp) for imgFp in glob.glob(imgDir+"/*.*")]
    mskFps = [os.path.join(mskDir, mskFp) for mskFp in glob.glob(mskDir+"/*.*")]
    assert len(imgFps)==len(mskFps)
    for i in range(len(imgFps)): f.write(f"{imgFps[i]},{mskFps[i]}\n") 
    print(f"Successfully generated csv at {saveFp}")

####################################
########### Computation ############
####################################

def computeDSC(pred, msk):
    return (2.*np.count_nonzero(np.multiply(pred, msk))+1)/(\
        np.count_nonzero(pred)+np.count_nonzero(msk)+1)

def computeBalanceWeights(mskDir, numClasses):
    """ Calculate the balanced weights based 
    on the number of voxels of each class 
    """
    mskFps = glob.glob(mskDir+"/*.*")
    weights = [0 for _ in range(numClasses)]
    print("Calculating the balance weights.....")
    for mskFp in mskFps:
        msk = preprocessMsk(readImg(mskFp))
        unique, counts = np.unique(msk, return_counts=True)
        weights = [weights[num]+counts[num] for num in range(numClasses)]
    assert 0 not in weights
    weights = [1./(weight) for weight in weights]
    print(f"Balanced factors: {weights}")
    return weights

####################################
############# 3D plot ##############
####################################

def plot_3d(image, threshold=-300):
    import matplotlib.pyplot as plt
    from skimage import measure, morphology
    from mpl_toolkits.mplot3d.art3d import Poly3DCollection
    p = image.transpose(2,1,0)
    p = p[:,:,::-1]
    verts, faces = measure.marching_cubes(p, threshold)
    fig = plt.figure(figsize=(10, 10))
    ax = fig.add_subplot(111, projection='3d')
    mesh = Poly3DCollection(verts[faces], alpha=0.1)
    face_color = [0.5, 0.5, 1]
    mesh.set_facecolor(face_color)
    ax.add_collection3d(mesh)
    ax.set_xlim(0, p.shape[0])
    ax.set_ylim(0, p.shape[1])
    ax.set_zlim(0, p.shape[2])
    plt.show()

####################################
########### Create GIF #############
####################################

def applyMsk(img, predMsk, gtMsk):
    """ Apply predMsk and gtMsk on the original img,
    Returns a 3D img with each pixel of 3 channels"""
    assert img.shape==predMsk.shape==gtMsk.shape, "Unmatched dimension"
    R, G, B = img.copy(), img.copy(), img.copy()
    R[np.where((predMsk==1)&(gtMsk!=1))] = 255
    R[np.where((predMsk!=1)&(gtMsk==1))] = 0
    R[np.where((predMsk==1)&(gtMsk==1))] = 255
    G[np.where((predMsk==1)&(gtMsk!=1))] = 0
    G[np.where((predMsk!=1)&(gtMsk==1))] = 255
    G[np.where((predMsk==1)&(gtMsk==1))] = 255
    B[np.where((predMsk==1)&(gtMsk!=1))] = 0
    B[np.where((predMsk!=1)&(gtMsk==1))] = 0
    B[np.where((predMsk==1)&(gtMsk==1))] = 0
    img = np.stack((R, G, B),0)
    img = np.moveaxis(img, 3, 0)
    return img


def applyContour(img, predMsk, gtMsk):
    assert img.shape==predMsk.shape==gtMsk.shape, "Unmatched dimension"
    img, R, G, B = img.copy(), img.copy(), img.copy(), img.copy()
    predMsk = predMsk.astype("uint8")
    gtMsk = gtMsk.astype("uint8")
    slices = []
    print("Computing contour array.....")
    for i in range(img.shape[2]):
        imgSlice = img[:,:,i]
        rgbSlice = np.stack((imgSlice,)*3,-1)
        predSlice  = predMsk[:,:,i]
        gtSlice  = gtMsk[:,:,i]
        _, predCts, _ = cv2.findContours(predSlice.copy(),\
            cv2.RETR_LIST,cv2.CHAIN_APPROX_SIMPLE)
        _, gtCts, _ = cv2.findContours(gtSlice.copy(),\
            cv2.RETR_LIST,cv2.CHAIN_APPROX_SIMPLE)
        if len(predCts)!=0:
            predCt = predCts[0]
            cv2.drawContours(rgbSlice, [predCt], 0, (255,0,0), 1)
        if len(gtCts)!=0:
            gtCt = gtCts[0]
            cv2.drawContours(rgbSlice, [gtCt], 0, (0,255,0), 1)
        slices.append(rgbSlice)
    img = np.asarray(slices)
    img = np.moveaxis(img, 3, 1)
    print("Generated contour array!")
    return img


def generateGIF(img, saveFp, fps):
    write_gif(img, saveFp, fps)
    print(f"Successfully wrote gif at: {saveFp}")


def isotropicDir(srcDir, destDir, zoomOrder):
    imgFps = glob.glob(srcDir+"/*.*")
    for imgFp in imgFps:
        saveFp = destDir + "/" + get_id(imgFp) + ".nii"
        isotropic(imgFp, saveFp=saveFp, zoomOrder=zoomOrder)


def getContourFromMsk(msk):
    newMsk = np.zeros(msk.shape)
    for i in range(msk.shape[2]):
        contours, _ = cv2.findContours(msk[:, :, i], cv2.RETR_LIST, cv2.CHAIN_APPROX_SIMPLE)
        newSlice = np.stack((newMsk[:, :, i],)*3, -1)
        cv2.drawContours(newSlice, contours, -1, (0, 255, 0), 1)
        newSlice = newSlice[:,:,1]
        newSlice[np.where(newSlice!=0)] = 1
        newMsk[:,:,i] = newSlice
    return newMsk