% Alberto Cereser, 2 Oct 2015
% Department of Physics, Technical University of Denmark
% alcer@fysik.dtu.dk

% Script to clean, for each projection, the set of combined images before
% the HT. Output:
% - Combined images, in Comb_blobs_final;
% - Datafile with the properties of the combined images (Data_blobs_final.txt);
% - Datafile keeping track of how the different images are combined
% (Image_combination_before_HT.txt)

clear; close all;

% Load data relative to the combined extinction spots, for all projections.
% File made using Extract_data_Comb_blobs.m directly from the images
Data = load('Data_blobs.txt');
rad = 10;    % Radius of the tube used to group nearby points

% Text file where we store the image name, so we know how they are combined
fid = fopen('Image_combination_before_HT.txt', 'wt');

for Omega = 1:6:181  % This is the projection number. The first proj
                    % corresponds to 0 degrees
    9_Funct_im_comb_2(Omega - 1, Data, rad, fid);
end
fclose(fid);

% Read in the data, directly from the combined images
imagefiles = dir(sprintf('Comb_blobs_final/'));
nfiles = length(imagefiles);    % Number of files found
path_dir = sprintf('Comb_blobs_final/');

Data_blobs = zeros(nfiles-3, 5);
for aa = 1:nfiles-3
    currentfilename = imagefiles(aa+3).name;
    path_image = strcat(path_dir, currentfilename);
    % Get image characteristics from the file name
    Im_data = sscanf(currentfilename,'%*4c_%*5c_%u_%u.png', 2);
    Omega = Im_data(1);
    IM_id = Im_data(2);
    A = imfill(double(imread(path_image)), 'holes');
    IM = int8(A);
    IM_scaled = IM./max(max(IM));
    Area = nnz(IM_scaled);
    centr = regionprops(IM_scaled,'centroid');
    centroids = cat(1, centr.Centroid);
    X_CM = centroids(1);
    Y_CM = centroids(2);
    Data_blobs(aa, :) = [Omega, IM_id, X_CM, Y_CM, Area];
end

% Save data relative to all the extinction spots
dlmwrite('Data_blobs_final.txt', Data_blobs, 'delimiter',' ', 'precision', 6);
