% Alberto Cereser, August 2015
% Technical University of Denmark

% This script locates the different extinction spots in a bw image, selects
% the ones with area bigger then a threshold value and saves them.

function 2_Blobs_detector_function(Omega, Lambda)
frame_size = 30;            % Size of the frame used to clean the image from noisy points
area_threshold = 100;       % We are interested in blobs with at least this number of pixels
area_threshold_max = 2000;  % Max num of pixels

filename = sprintf('../png/Im_div_roll_median_Fe_%03i_%05i_filt_bin.png', Omega-1, Lambda);
IM = double((imread(filename))');
[x_IM, y_IM] = size(IM);

% Make a folder where to store the isolated blobs
fn = fullfile(sprintf('Isolated_blobs_%03i', area_threshold));
if exist(fn, 'dir')
    warning('The output folder exists');
else
    mkdir(fn);
end

% To avoid noise, we clear the image using a frame
IM_clean = zeros(x_IM, y_IM);
for i = frame_size:x_IM-frame_size
    for j = frame_size:y_IM-frame_size
        IM_clean(i,j) = IM(i,j);
    end
end

[L,num_blobs] = bwlabel(IM_clean);  % Label the blobs in the cleaned image

% Select the big blobs (with area > area_threshold) and plot them.
Area = zeros(num_blobs, 1);
n_big_blobs = 0;
for k = 1:num_blobs
    Area(k) = sum(sum(L==k));
    if Area(k) > area_threshold && Area(k) < area_threshold_max
        n_big_blobs = n_big_blobs + 1;
    end
end

list_big = zeros(n_big_blobs, 2);
l = 0;
for k = 1:num_blobs
    if Area(k) > area_threshold && Area(k) < area_threshold_max
        l = l + 1;
        list_big(l,1) = k;
        list_big(l,2) = Area(k);
    end
end

disp('The big blobs have these tags:')
disp((list_big(:,1))');

L_big = zeros(x_IM, y_IM);
for m = 1:n_big_blobs               % To select a blob, select a value for m
    for ii = 1:x_IM
        for jj = 1:y_IM
            if L(ii,jj) == list_big(m)
                L_big(ii,jj) = L(ii,jj);
            end
        end
    end
end

% Separately save an image for each big blob
for m = 1:n_big_blobs
    tag = list_big(m,1);
    Im_tag = zeros(x_IM, y_IM);
    for ii = 1:x_IM
        for jj = 1:y_IM
            if L(ii,jj) == list_big(m)
                Im_tag(ii,jj) = 1;
            end
        end
    end
    Im_tag = imfill(Im_tag);
    output = sprintf('Isolated_blobs_%03i/Isolated_blob_%03i_%05i_%02i.png', area_threshold, Omega-1, Lambda, m);
    imwrite(Im_tag', output);
end

% Return a table with the characteristics of the big blobs
T = regionprops('table',L_big,'Centroid',...
    'MajorAxisLength','MinorAxisLength', 'Area');   % More options at http://se.mathworks.com/help/images/ref/regionprops.html

TF = ismissing(T);
T2 = T(~any(TF,2),:)
