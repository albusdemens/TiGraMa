% Alberto Cereser, 13 October 2015
% Technical University of Denmark, alcer@fysik.dtu.dk

% This script extracts, from the blobs in Comb_blobs/, combined using
% 5_Track_lambda.m, a file with their main characteristics (Omega, ID, coord
% centroid, area)
% To divide in separate folders images relative to different proj, use
% 6_commands_sim_diff_lambda.sh

% Let's start by reading in the data, directly from the combined images
imagefiles = dir(sprintf('/Comb_blobs/'));
nfiles = length(imagefiles);    % Number of files found
path_dir = sprintf('/Comb_blobs/');

Data_blobs = zeros(nfiles-3, 5);
for aa = 1:nfiles-3
    currentfilename = imagefiles(aa+3).name;
    path_image = strcat(path_dir, currentfilename);
    % Get image characteristics from the file name
    Im_data = sscanf(currentfilename,'%*4c_%*4c_%u_%u.png', 2);
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
dlmwrite('Data_blobs.txt', Data_blobs, 'delimiter',' ', 'precision', 6);
