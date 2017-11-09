% Alberto Cereser, 20 October 2015
% Technical University of Denmark, alcer@fysik.dtu.dk

% Script to plot in 3D data the slices obtained from Voxels_tagging_final.m
% Just an alternative to the standard Fiji procedure

clear; close all;

% Let's start by reading in the data, directly from the sample layers
imagefiles = dir(sprintf('Grain_reconstruction/'));
nfiles = length(imagefiles);    % Number of layers
path_dir = sprintf('Grain_reconstruction/');
Data_blobs = zeros(nfiles-3, 5);
Voxellized_sample = zeros(200,200,nfiles);
Voxellized_sample_edit = zeros(200,200,nfiles);

for aa = 1:nfiles-3
    currentfilename = imagefiles(aa+3).name;
    path_image = strcat(path_dir, currentfilename);
    A = double(imread(path_image));
    Voxellized_sample(:,:,aa) = A(:,:);
end

% To plot the sample with PATCH_3Darray, we need the empty voxels to have
% value nan
for ii = 1:200
    for jj = 1:200
        for kk = 1:nfiles
            if Voxellized_sample(ii,jj,kk) == 0
                Voxellized_sample_edit(ii,jj,kk) = NaN;
            else
                Voxellized_sample_edit(ii,jj,kk) = Voxellized_sample(ii,jj,kk);
            end
        end
    end
end

hpat = PATCH_3Darray(Voxellized_sample_edit, 'col');

% Rotate sample and save result as a video
OptionZ.FrameRate=15;OptionZ.Duration=5.5;OptionZ.Periodic=true;
CaptureFigVid([-20,10;-110,10;-190,80;-290,10;-380,10], 'WellMadeVid',OptionZ)
