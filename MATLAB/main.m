% Han - generate mesh
clc;
clear all;
close all;

% input arguments
niiImgFp = "/home/liuh26/Desktop/trajactory-data/thalamus/Normal003-right.nii";
saveDir = "/home/liuh26/Desktop";
thresh = 0.5;


% convert nifti to mesh
nii2mesh(niiImgFp, thresh, saveDir);
