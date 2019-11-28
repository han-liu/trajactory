function obj = nii2surfstruct(imgFp, thresh, saveFlag)
% Convert an Nifti image to a surfstruct object

% Input arguments:
% imgFp: file path of the input image
% thresh: the threshold used to generate the isosurface
% saveFlag: logical value. i.e. true/false. Determine if to save the output
%           into mat file under the same directory of the input image.

% Output:
% obj: the output surfstruct object

% Assumption:
% 1. Isotropic
% 2. Orientation is [1;2;3]

V = niftiread(imgFp);
[m,n,p] = size(V);
[X,Y,Z] = meshgrid(1:n,1:m,1:p);
[triangles, vertices] = isosurface(X,Y,Z, V, thresh);
dim = [m;n;p];
vox = [1;1;1];
orient = [1;2;3];
color = [255;255;255];
obj = surfstruct(vertices, triangles, dim, vox, orient, color);
if saveFlag
    save('surfstruct_obj.mat', 'obj')
end
fprintf("** Successfully saved the surfstruct object **\n")
end