% Alberto Cereser, 18 Nov 2015
% Department of Physics, Technical University of Denmark
% alcer@fysik.dtu.dk

% Function to
% a) find the peaks in the accumulator matrix;
% b) plot the corresponding fitting function

function [R_fit, alpha_fit] = 15_Funct_fit_HT(Counter_matrix1, threshold, ...
    n_large_pixels, L_max, OL, threshold_distance_fit)

Number_points = size(OL, 1);
% Find center of the curves in the Hough transform. To do so, I binarize
% the image using the max
[dim_X, dim_Y] = size(Counter_matrix1);
hist_bin = zeros(dim_X, dim_Y);
hist_bin_clean = zeros(dim_X, dim_Y);
max_hist1 = max(max(Counter_matrix1));
for qq = 1:dim_X
    for rr = 1:dim_Y
        if Counter_matrix1(qq,rr) > max_hist1*0.6
            hist_bin(qq,rr) = Counter_matrix1(qq,rr);
        end
    end
end

% Image analysis: erosion, dilation (excluded after missing one maxima),
% take out point outside, find CM
se = strel('disk',1);
Image_final = imclose(hist_bin,se);
[L,num] = bwlabel(Image_final);
% Let's exclude the single points
counts = sum(bsxfun(@eq,L(:),1:num));
valid_blobs = find(counts>threshold);
number_valid = length(find(counts>threshold));

if (number_valid~=0)
    for blob_num = (1:number_valid)
        for x = (1:dim_X)
            for y = (1:dim_X)
                if (L(x,y) == valid_blobs(blob_num))
                    hist_bin_clean(x,y) = 1;%Image_final(x,y);
                end
            end
        end
    end
end

% The Matlab function to find peaks in the Hough transform may find
% more peaks than expected; better to use my own function

Ilabel = bwlabel(hist_bin_clean);
stat = regionprops(Ilabel,'Centroid');

number_centers = numel(stat);
coord_centers = zeros(number_centers,2);
for x = 1:number_centers
    coord_centers(x,1) = stat(x).Centroid(1);   % Alpha
    coord_centers(x,2) = stat(x).Centroid(2);   % Radius
end

for xx = 1:number_centers
    alpha_fit(xx) = coord_centers(xx,1)*360/n_large_pixels;
    R_fit(xx) = coord_centers(xx,2)*L_max/n_large_pixels;
    % Depending on the binning value, the sin function(s) could be
    % translated. Let's check that
    % alpha_fit1(xx) = coord_centers(xx,1)*360/n_large_pixels-10;
    % R_fit1(xx) = coord_centers(xx,2)*5/n_large_pixels;
end

if number_centers == 0
    alpha_fit(1) = 0;
    R_fit(1) = 0;
end

% For each point, we want to find what is the closest fitting line.
% This is done in two steps: at first we check if the distance is below
% a threshold, and then we look for the lowest distance value.
fit_value = zeros(Number_points, number_centers);
dist_curve = zeros(Number_points, number_centers);
dist_curve_thresh = zeros(Number_points, number_centers);
% Here we store the data after using a threshold
dist_curve_thresh_fin = zeros(Number_points, number_centers);
% Here we store the data after calculating the minimal distance
number_candidates = zeros(Number_points);
for aa = 1:Number_points
    for bb = 1:number_centers
        fit_value(aa,bb) = double(R_fit(bb)*cosd(alpha_fit(bb) + OL(aa,2)));
        dist_curve(aa,bb) = abs(fit_value(aa,bb) - OL(aa,4));
        if (dist_curve(aa,bb) < threshold_distance_fit)
            dist_curve_thresh(aa,bb) = dist_curve(aa,bb);
        end
    end
end
for cc = 1:Number_points
    for dd = 1:number_centers
        number_candidates(cc) = sum(dist_curve_thresh(cc,:)~=0,2);
        if number_candidates(cc) > 1
            index_min_dist = find(dist_curve_thresh(cc, :) > 0, 1, 'first');
            dist_curve_thresh_fin(cc,index_min_dist) = dist_curve_thresh(cc, index_min_dist);
        else
            dist_curve_thresh_fin(cc,dd) = dist_curve_thresh(cc,dd);
        end
    end
end

% We label each point in the Angle, Y _measured space so that we know
% what's the closest curve
tagged_points = zeros(Number_points,3);
tagged_points_edit = zeros(Number_points, number_centers + 1);
for ee = 1:Number_points
    for ff = 1:number_centers
        if dist_curve_thresh_fin(ee,ff) > 0
            tagged_points(ee,1) =  OL(ee,2);
            tagged_points(ee,2) =  OL(ee,4);
            tagged_points(ee,3) =  ff;
        end
    end
    for gg = 1:number_centers
        if tagged_points(ee,3) == gg
            tagged_points_edit(ee,1) = tagged_points(ee,1); % Angle
            tagged_points_edit(ee,gg+1) = tagged_points(ee,2);  % Y_measured
        end
    end
end

% To optimize the threshold value, let's check how many empty rows we
% have in dist_curve_thresh_fin
percentage_labelled = nnz(dist_curve_thresh_fin)/Number_points;
