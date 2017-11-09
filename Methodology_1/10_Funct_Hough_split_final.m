% Alberto Cereser, 6 Jun 2015
% Technical University of Denmark, alcer@fysik.dtu.dk

% Function to calculate the Hough transform for a given set of points. Function
% called in Hough_trans_pol_section_Murofi_split.m

function [number_centers, tagged_points_edit, alpha_fit, R_fit] = ...
    10_Funct_Hough_split_final(Y_measured, Angle, center, ...
    n_large_pixels, Number_points, steps_alpha, steps_size, threshold, ...
    threshold_distance_fit)

Y_measured_correct = Y_measured - center; %With this correction, we take into
% account the position of the cylinder's axis, calculated using
% Hough_transf_RC.m

Counter_matrix = zeros(n_large_pixels, n_large_pixels); % We use this matrix to locate the maxima
% of the Hough transform
Counter_matrix_bin = zeros(n_large_pixels, n_large_pixels); % This is the binarized version of
% Counter_matrix, made using the threshold "cut"
Counter_matrix_one_blob_sum = zeros(n_large_pixels, n_large_pixels);

for pp = (1:Number_points)
    Counter_matrix_one_blob = zeros(n_large_pixels, n_large_pixels);
    R_n = zeros(1,(steps_alpha/steps_size));
    A_n = zeros(1,(steps_alpha/steps_size));
    n = 0;
    for Alpha = (0:steps_size:steps_alpha-1)  % This is the angle of the blob CM in the sample
        den = sind(Angle(pp)+Alpha);
        if ((abs(den)>2*abs(Y_measured_correct(pp))*0.055/10) && (Y_measured_correct(pp)*0.055/den > 0))  % Conditions included to avoid problems with the den
            % In the formula above, L is the diameter of the cylinder
            % (L; value explicited to avoid bug)
            R = (Y_measured_correct(pp))*0.055/den;
            n = n + 1;
            R_n(n) = R;
            A_n(n) = Alpha;
        end
    end
    len_n_one_alpha = size(A_n,2);
    M_n_one_alpha = zeros(n_large_pixels,n_large_pixels);
    for ij = 1:n
        ix = floor(A_n(ij)/360*n_large_pixels)+1;
        iy = floor(R_n(ij)/5*n_large_pixels)+1;
        M_n_one_alpha(ix,iy) = 1;
    end

    for ii = 1:n_large_pixels
        for jj = 1:n_large_pixels
            Counter_matrix_one_blob_sum(ii,jj) = Counter_matrix_one_blob_sum(ii,jj) + M_n_one_alpha(ii,jj);
        end
    end
end

Counter_matrix1 = Counter_matrix_one_blob_sum';
%Counter_matrix1(size(n,1) + 1, size(n,2) + 1) = 0;

%g = figure();
%h = pcolor(Counter_matrix1); title('Counter matrix');
%g1 =  figure();
%h1 = mesh(Counter_matrix1);

% Find center of the curves in the Hough transform. To do so, I binarize
% the image using the max
[dim_X, dim_Y] = size(Counter_matrix1);
hist_bin = zeros(dim_X, dim_Y);
hist_bin_clean = zeros(dim_X, dim_Y);
max_hist1 = max(max(Counter_matrix1));
for qq = 1:dim_X
    for rr = 1:dim_Y
        if Counter_matrix1(qq,rr) > max_hist1*0.5
            hist_bin(qq,rr) = Counter_matrix1(qq,rr);
        end
    end
end

% Image analysis: erosion, dilation (excluded after missing one maxima),
% take out point outside, find CM
se = strel('disk',3);
Image_final = imclose(hist_bin,se);
[L,num] = bwlabel(Image_final);
% Let's exclude the single points
counts = sum(bsxfun(@eq,L(:),1:num));
valid_blobs = find(counts>threshold);
number_valid = length(find(counts>threshold));

%g2 = figure();
%h2 = mesh(L);

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

%g2 = figure();
%h2 = mesh(hist_bin_clean);

% The Matlab function to find peaks in the Hough transform may find
% more peaks than expected; better to use my own function

Ilabel = bwlabel(hist_bin_clean);
stat = regionprops(Ilabel,'Centroid');

number_centers = numel(stat);
coord_centers = zeros(number_centers,2);
xy_cylinder = zeros(number_centers,2);
for x = 1:number_centers
    coord_centers(x,1) = stat(x).Centroid(1);   % Alpha
    coord_centers(x,2) = stat(x).Centroid(2);   % Radius
    xy_cylinder(x,1) = coord_centers(x,2)*cosd(coord_centers(x,1));
    xy_cylinder(x,2) = coord_centers(x,2)*sind(coord_centers(x,1));
end

for xx = 1:number_centers
    alpha_fit(xx) = coord_centers(xx,1)*360/n_large_pixels;
    R_fit(xx) = coord_centers(xx,2)*5/n_large_pixels;
    % Depending on the binning value, the sin function(s) could be
    % translated. Let's check that
    alpha_fit1(xx) = coord_centers(xx,1)*360/n_large_pixels-10;
    R_fit1(xx) = coord_centers(xx,2)*5/n_large_pixels;
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
        angular_coord = Angle(aa);
        fit_value(aa,bb) = double((R_fit(bb)*sind(alpha_fit(bb)+angular_coord)/0.055) + 48);
        dist_curve(aa,bb) = abs(fit_value(aa,bb) - Y_measured(aa));
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
            tagged_points(ee,1) =  Angle(ee);
            tagged_points(ee,2) =  Y_measured(ee);
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
