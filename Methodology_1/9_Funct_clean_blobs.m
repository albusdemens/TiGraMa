% Alberto Cereser, 2 October 2015
% Technical University of Denmark, alcer@fysik.dtu.dk

% Simple function that fills holes in a blob and cleans disconnected
% pixels. To be used after combining images

function [biggest_blob, Area_comb, X_CM, Y_CM] = Funct_clean_blobs(IM)

% Fill holes and avoid disconnected pixels
IM_filled = imfill(IM, 'holes');
[L,num] = bwlabel(IM_filled);
count_pixels_per_obj = sum(bsxfun(@eq,L(:),1:num));
[~,ind] = max(count_pixels_per_obj)

if isempty(ind) == 0
    biggest_blob = double(L==ind);

    Area_comb = nnz(biggest_blob);
    centr = regionprops(biggest_blob,'centroid');
    centroids = cat(1, centr.Centroid);
    X_CM = centroids(1);
    Y_CM = centroids(2);
else
    biggest_blob = 0;
    Area_comb = 0;
    X_CM = 0;
    Y_CM = 0;
end
