import utils
import math
import os
import numpy as np
from deepbrain import Extractor
import scipy.ndimage.morphology



# generate masks of different tissues from SLANT output
def generateSegFromSlant(slantMskFp, saveDir):
    # SLANT segmentation labels:
    # 59: Right thalamus
    # 60: Left thalamus
    # 4: 3rd ventricle
    # 11: 4th ventricle
    # 51: Right Lateral Ventricle    # 52: Left Lateral Ventricle
    names = ["right-thalamus", "left-thalamus", "3rd-ventricle", "4th-ventricle", "RLV", "LLV"]
    for i, label in enumerate([59, 60, 4, 11, 51, 52]):
        msk = utils.readImg(slantMskFp)
        msk[np.where(msk != label)] = 0
        msk[np.where(msk != 0)] = 1
        msk = msk.astype(dtype="uint16")
        utils.saveImg(msk, saveDir+"/Normal003-"+names[i]+".nii")


# generate brain-surface
def generateBrainSurf(rawImgFp, saveDir, thresh=15):
    img = utils.readImg(rawImgFp)
    img[np.where(img < thresh)] = 0  # global threshold
    img = scipy.ndimage.morphology.binary_erosion(img)  # erosion
    img = utils.getLargestConnectedComp(img)  # largest connected component
    img = scipy.ndimage.morphology.binary_dilation(img)  # dilation
    # eliminate the outlier
    for y in range(14):
        for z in range(114, 119):
            for x in range(256):
                img[x, y, z] = 0
    img = img.astype(dtype="uint16")
    utils.saveImg(img, saveDir + "/brain-surf.nii")
    return img


# generate target brain-surface
def generateTargetSurf(brainSurfImgFp, saveDir):
    img = utils.readImg(brainSurfImgFp)
    for x in range(256):  # Anterior-Posterior
        for y in range(176):  # Left-Right
            for z in range(176):  # Superior-Inferior
                # define the boundaries of x, y and z
                if x < 60 or x > 117 or z < 100:
                    img[x, y, z] = 0
    img = img.astype(dtype="uint16")
    utils.saveImg(img, saveDir + "/target-surf.nii")
    return img


# fill the surface
def fillSurf(targetSurfImgFp, saveDir):
    img = utils.readImg(targetSurfImgFp)
    for y in range(176):
        for x in range(60, 117):
            top = 0
            for z in range(100, 176)[::-1]:
                if img[x, y, z] == 1:
                    top = z
                    break
            for z in range(100, top):
                img[x, y, z] = 1
    img = scipy.ndimage.morphology.binary_dilation(img)
    img = img.astype(dtype="uint16")
    utils.saveImg(img, saveDir + "/fill-surf.nii")
    return img


# generate the surface
def generateSurf(targetSurfImgFp, saveDir):
    img = utils.readImg(targetSurfImgFp)
    img = scipy.ndimage.morphology.binary_dilation(img)
    surf = np.zeros(img.shape)
    for x in range(60, 117):
        for y in range(176):
            # from superior to inferior
            for z in range(100,176)[::-1]:
                if img[x, y, z] == 1:
                    surf[x, y, z] = 1
                    break
    for x in range(60, 117):
        for z in range(100, 176):
            # from right to left
            for y in range(176)[::-1]:
                if img[x, y, z] == 1:
                    surf[x, y, z] = 1
                    break
            # from left to right
            for y in range(176):
                if img[x, y, z] == 1:
                    surf[x, y, z] = 1
                    break
    # for y in range(176):
    #     for z in range(100, 176):
    #         # from anterior to posterior
    #         for x in range(60, 117):
    #             if img[x, y, z] == 1:
    #                 surf[x, y, z] = 1
    #                 break
    surf = surf.astype(dtype="uint16")
    utils.saveImg(surf, saveDir+"/surface.nii")
    return surf


# post-processing surface
def postProcessSurf(surfImgFp, saveDir):
    img = utils.readImg(surfImgFp)
    img = scipy.ndimage.morphology.binary_dilation(img)
    img = img.astype(dtype="uint16")
    utils.saveImg(img, saveDir + "/post-surface.nii")
    return img


# skull-stripped. Returns a binary mask
def skullStrip(imgFp, saveDir):
    img = utils.readImg(imgFp)
    ext = Extractor()
    prob = ext.run(img)
    mask = (prob > 0.5)
    mask = mask.astype(dtype="uint16")
    utils.saveImg(mask, saveDir + "/skull-strip.nii")
    return mask


# visualize the trajectory
# EXAMPLE:
#           startingPt = np.array([164, 87, 63])
#           targetPt   = np.array([93, 123, 99])
#  In C++ format, but in y, x ,z order in Python

def visualizeTraj(postSurfImgFp, vesselImgFp, ventricleImgFp, brainSurfImgFp, ANTImgFp, startingPt, targetPt, saveDir):
    postSurf = utils.readImg(postSurfImgFp)
    dist = ((targetPt[0] - startingPt[0]) ** 2 + (targetPt[1] - startingPt[1]) ** 2
            + (targetPt[2] - startingPt[2]) ** 2) ** 0.5
    flow = [(targetPt[0] - startingPt[0]) / dist, (targetPt[1] - startingPt[1]) / dist,
            (targetPt[2] - startingPt[2]) / dist]
    traj = np.zeros(postSurf.shape)
    nPt = int(dist)
    print(f"Direction of the trajectory: ", flow)
    print(f"Number of points: ", nPt)
    print("** Generating the trajectory **")
    for idx in range(nPt):
        x = int(round(startingPt[1] + flow[1] * idx))
        y = int(round(startingPt[0] + flow[0] * idx))
        z = int(round(startingPt[2] + flow[2] * idx))
        traj[x, y, z] = 1

    # extend the trajectory
    for idx in range(25):
        x = int(round(startingPt[1] - flow[1] * idx))
        y = int(round(startingPt[0] - flow[0] * idx))
        z = int(round(startingPt[2] - flow[2] * idx))
        traj[x, y, z] = 1

    traj = scipy.ndimage.morphology.binary_dilation(traj)
    traj = traj.astype(dtype="uint16")

    # plot the target point
    target = np.zeros(postSurf.shape)
    target[targetPt[1]-2:targetPt[1]+2, targetPt[0]-2:targetPt[0]+2, targetPt[2]-2:targetPt[2]+2] = 1

    # include other components in our window
    vessel = utils.readImg(vesselImgFp)
    ventricle = utils.readImg(ventricleImgFp)
    brain = utils.readImg(brainSurfImgFp)
    ANT = utils.readImg(ANTImgFp)

    # Assign with different colors
    traj[np.where(traj != 0)] = 2000          # trajectory
    ventricle[np.where(ventricle != 0)] = 20  # ventricles
    postSurf[np.where(postSurf != 0)] = 10    # post-surface
    brain[np.where(brain != 0)] = 5           # brain-surface
    vessel[np.where(vessel != 0)] = 200       # vessel
    target[np.where(target != 0)] = 50        # target point
    ANT[np.where(ANT != 0)] = 300             # ANT

    # plot the final planning graph
    finalPlanning = ventricle + postSurf + traj + vessel + brain + ANT  # + target
    finalPlanning[np.where(finalPlanning == 15)] = 10
    finalPlanning[np.where(finalPlanning == 210)] = 10
    finalPlanning[np.where(finalPlanning == 215)] = 10
    finalPlanning[np.where(finalPlanning == 2000)] = 2005
    finalPlanning[np.where(finalPlanning == 2015)] = 2005
    finalPlanning = finalPlanning.astype(dtype="uint16")
    utils.saveImg(finalPlanning, saveDir+"/final_planning.nii")
    return finalPlanning


if __name__ == "__main__":
    # The entire pipeline of visualization

    # define the directories
    saveDir = "/home/liuh26/Desktop/new"
    if not os.path.exists(saveDir):
        os.mkdir(saveDir)

    # Pre-defined paths
    vesselImgFp = "/home/liuh26/Desktop/vessel.nii"
    ventricleImgFp = "/home/liuh26/Desktop/ventricles.nii"
    rawImgFp = "/home/liuh26/Desktop/Normal003.nii"
    slantMskFp = "/home/liuh26/Desktop/trajactory-data/SLANT/Normal003_seg.nii"
    ANTImgFp = "/home/liuh26/Desktop/ANT.nii"

    # To be generated in saving directory
    brainSurfImgFp = saveDir + "/brain-surf.nii"
    targetSurfImgFp = saveDir + "/target-surf.nii"
    fillSurfImgFp = saveDir + "/fill-surf.nii"
    surfImgFp = saveDir + "/surface.nii"
    postSurfImgFp = saveDir + "/post-surface.nii"

    # import the result from the C++ program
    # startingPt = np.array([164, 87, 63])
    startingPt = np.array([151, 114, 128])  # new vessel map
    targetPt   = np.array([92, 122, 98])

    # running the pipeline
    # generateSegFromSlant(slantMskFp, saveDir)
    generateBrainSurf(rawImgFp, saveDir)
    generateTargetSurf(brainSurfImgFp, saveDir)
    # fillSurf(targetSurfImgFp, saveDir)
    generateSurf(targetSurfImgFp, saveDir)
    postProcessSurf(surfImgFp, saveDir)
    visualizeTraj(postSurfImgFp, vesselImgFp, ventricleImgFp, brainSurfImgFp, ANTImgFp, startingPt, targetPt, saveDir)