% Alberto Cereser, 8 Jun 2015
% Department of Physics, Technical University of Denmark
% alcer@fysik.dtu.dk

% Developed for Hough_trans_pol_rolling_Murofi.m

function Function_HT_pol_rolling_final(vv, step_size, Data_blobs)
num_el = numel(Data_blobs(:,1));
center = 48;

output = 'CM_alpha_R.txt';
fid = fopen(output, 'a+');

%%%%%% PARAMETERS USED %%%%%%
threshold = 5;
% This is the minimum number of pixels to consider to calculate the
% position of the peak centers in the Hough transform
threshold_distance_fit = 5;    % VALUE CHANGED; original = 10
% Using this distance from the fit line, we group the points around the
% calculated curves
cut = 6;
% We use this value to cut the peak values in Counter_matrix
min_points = 7;
% To choose which results of the Hough transform to store, we check that at
% least this number of values is in a threshold_distance_fit band from
% the fitting line

IM_id = Data_blobs(:,2);
X = Data_blobs(:,4); % BL18 data 2015
Y = Data_blobs(:,3);
x_min = min(X);
x_max = max(X);
y_min = min(Y);
y_max = max(Y);

interval_num = 0;

Tagged_matrix = zeros(num_el, 6);
interval_num = interval_num + 1;
for zz = 1:num_el
    if ((X(zz) > (vv - (step_size-1)/2)) && (X(zz) < (vv + (step_size-1)/2)))
        Tagged_matrix(zz, 1) = X(zz);
        Tagged_matrix(zz, 2) = Y(zz);
        Tagged_matrix(zz, 3) = Data_blobs(zz, 1);   % projection
        Tagged_matrix(zz, 4) = interval_num;        % interval tag
        Tagged_matrix(zz, 5) = Data_blobs(zz, 5);   % Area
        Tagged_matrix(zz, 6) = IM_id(zz);
    end
end
tagged_lines = nnz(Tagged_matrix(1, :));
Tagged_matrix_clean = zeros(tagged_lines, 6);
clean_count = 0;
for uu = 1:num_el
    if Tagged_matrix(uu, 5) > 0
        clean_count = clean_count + 1;
        Tagged_matrix_clean(clean_count, :) = Tagged_matrix(uu, :);
        % Y, X, projection, tag, Area, lambda
    end
end
disp(vv);
[Data_with_CM] = Funct_HT_sec_Murofi_final(...
    Tagged_matrix_clean,threshold,threshold_distance_fit,...
    min_points,y_min, y_max, vv);
dlmwrite(output, Data_with_CM,'-append','delimiter','\t','precision',5);
fclose(fid);
