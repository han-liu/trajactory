#include <vector>
#include <iostream>
#include <itkImage.h>
#include <itkImageFileWriter.h>
#include <itkImageFileReader.h>
#include <itkSignedMaurerDistanceMapImageFilter.h>
#include "itkLinearInterpolateImageFunction.h"
#include "itkContinuousIndex.h"
#include <cmath>
#include "itkRescaleIntensityImageFilter.h"


using ImageType = itk::Image<double, 3>;
using ReaderType = itk::ImageFileReader<ImageType>;
using WriterType = itk::ImageFileWriter<ImageType>;
using SignedMaurerDistanceMapImageFilterType = itk::SignedMaurerDistanceMapImageFilter<ImageType, ImageType>;
using RescaleType = itk::RescaleIntensityImageFilter<ImageType, ImageType>;
using IteratorType = itk::ImageRegionIterator<ImageType>;


// compute the minimum distance given a trajectory (defined by a starting point and a target point) and a distance map.
double computeMinDist(std::vector<int> startingPt, std::vector<int> targetPt, ImageType::Pointer distMap){
    double dist = sqrt(pow(startingPt[0]-targetPt[0],2)+pow(startingPt[1]-targetPt[1],2)+pow(startingPt[2]-targetPt[2],2));
    int nPt = int(dist);
    std::vector<double> flow = {0, 0, 0};
    for (int i=0; i<3; i++) {flow[i] = (targetPt[i] - startingPt[i])/double(nPt);}
    double minDist = 1e3;
    for (int j=0; j<nPt; j++){
        itk::LinearInterpolateImageFunction<ImageType, double>::Pointer interpolator = itk::LinearInterpolateImageFunction<ImageType, double>::New();
        itk::ContinuousIndex<double, 3> pixel;
        for (int u=0; u<3; u++) {pixel[u] = startingPt[u] + flow[u]*j;}
        interpolator->SetInputImage(distMap);
        double d = interpolator->EvaluateAtContinuousIndex(pixel);
        if ( d <= 0){return 0;}
        else { if (d < minDist){ minDist = d; }}
    }
    return minDist;
}


int main(int argc, char* argv[])
{
    std::cerr << "Usage:\n[1] binary surface image\n[2] ventricle mask\n[3] vessel mask\n[4] ventricle distance map\n[5] vessel distance map\n" << std::endl;

    // initialize
    ImageType::Pointer surf = ImageType::New();
    ImageType::Pointer venMask = ImageType::New();
    ImageType::Pointer vesMask = ImageType::New();
    ImageType::Pointer venDistMap = ImageType::New();
    ImageType::Pointer vesDistMap = ImageType::New();

    // read the input images
    ReaderType::Pointer surfReader = ReaderType::New();
    ReaderType::Pointer venReader = ReaderType::New();
    ReaderType::Pointer vesReader = ReaderType::New();

    surfReader->SetFileName(argv[1]);
    surfReader->Update();
    surf = surfReader->GetOutput();

    venReader->SetFileName(argv[2]);
    venReader->Update();
    venMask = venReader->GetOutput();

    vesReader->SetFileName(argv[3]);
    vesReader->Update();
    vesMask = vesReader->GetOutput();

    // generate distance map of ventricles
    std::cout << "** Generating distance map of ventricles **" << std::endl;
    SignedMaurerDistanceMapImageFilterType::Pointer venDistFilter = SignedMaurerDistanceMapImageFilterType::New();
    WriterType::Pointer venDistWriter = WriterType::New();
    venDistFilter->SetInput(venMask);
    venDistMap = venDistFilter->GetOutput();
    venDistWriter->SetFileName(argv[4]);
    venDistWriter->SetInput(venDistMap);
    venDistWriter->Update();

    // generate distance map of vessels
    std::cout << "** Generating distance map of vessels **" << std::endl;
    SignedMaurerDistanceMapImageFilterType::Pointer vesDistFilter = SignedMaurerDistanceMapImageFilterType::New();
    WriterType::Pointer vesDistWriter = WriterType::New();
    vesDistFilter->SetInput(vesMask);
    vesDistMap = vesDistFilter->GetOutput();
    vesDistWriter->SetFileName(argv[5]);
    vesDistWriter->SetInput(vesDistMap);
    vesDistWriter->Update();

    IteratorType surfIterator(surf, surf->GetLargestPossibleRegion());
    IteratorType venIterator(venMask, venMask->GetLargestPossibleRegion());
    IteratorType vesIterator(vesMask, vesMask->GetLargestPossibleRegion());
    std::vector<int> targetPt = {93, 123, 99};

    std::cout << "** Computing the optimal trajectory **" << std::endl;
    double opt_dist = 0;
    std::vector<int> opt_startingPt = {0, 0, 0};

    for (surfIterator.GoToBegin(), venIterator.GoToBegin(), vesIterator.GoToBegin(); !surfIterator.IsAtEnd(); ++venIterator, ++vesIterator, ++surfIterator){
        if (surfIterator.Get() == 1){
            ImageType::IndexType idx = surfIterator.GetIndex();
            std::vector<int> startingPt = {int(idx[0]), int(idx[1]), int(idx[2])};
//            std::cout << "starting Pt index: " << idx[0] << " " << idx[1] << " " << idx[2] << std::endl;
            double venDist = computeMinDist(startingPt, targetPt, venDistMap);
            double vesDist = computeMinDist(startingPt, targetPt, vesDistMap);
            if (venDist > opt_dist) {
                std::cout << "Found better trajectory at ";
                opt_dist = venDist;
                for (int i = 0; i < 3; i++) {opt_startingPt[i] = startingPt[i]; std::cout << startingPt[i] << " ";}
                std::cout << " minDist: " << venDist;
                std::cout << std::endl;
            }
        if (venDist+vesDist > opt_dist) {
            opt_dist = venDist + vesDist;
            for (int i = 0; i < 3; i++) { opt_startingPt[i] = startingPt[i]; }
        }
        }
    }

    std::cout << "** Finished **" << std::endl;
    return 0;
}
