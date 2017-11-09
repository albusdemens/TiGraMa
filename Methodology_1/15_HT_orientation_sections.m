% Alberto Cereser and Martina Zamboni, 16 Nov 2015
% Department of Physics, Technical University of Denmark
% alcer@fysik.dtu.dk

% Hough transform code used to calculate the grain orientation by finding curves
% in the (Omega, lambda) plot

clear; close all;

%step_size = 25;    % Should be an odd number
grain_num = 3;      % Grain number in the OL list

n_large_pixels = 75;   % Size of the mesh used for the HT. Lambda > 2.5 AA
n_large_pixels_2 = 100; % Lambda in 2, 3 AA
n_large_pixels_3 = 200; % Lambda in 1.6, 2.1 AA
steps_alpha = 360;
steps_size = 0.2; % 0.1

lambda_min = 2.5;
cut_proj = 0.8;     % The only projections considered are those with # el <
                    % max*cut
l_shift = 0.0731;   % Number from Study_lambda_shift.m, where the distance
                    % Expected peaks - exp values is calculated
threshold = 1;      % This is the minimum number of pixels to consider to
                    % calculate the position of the peak centers in the
                    % Hough transform. Lambda > 2.5 AA
threshold_2 = 2;    % Lambda in 1.5, 3 AA
threshold_3 = 3;
threshold_distance_fit = 0.15;
                    % Using this distance from the fit line, we group the
                    % points around the calculated curves
L_max = 4.0538;     % Maximum lambda value (110 Bragg edge)
L_shift = 0.0731;   % Lambda shift, calculated using Study_lambda_shift.m

% Start by reading in the data
OL = load('OL_final_OK.txt');
n_grains = max(OL(:,1));

[r, ] = find(OL(:,5) == grain_num);   % Values relative to a single grain
OL_g(:,1:3) = OL(r,1:3);
OL_g(:,4) = OL(r,4) + L_shift;

% Bin the data. Parameter found using Map_bin_parameters.m
Omega_val = unique(OL_g(:,2));
n_bin = 0;
for hh = 1:size(Omega_val,1)
    [r_o,] = find(OL_g(:,2) == Omega_val(hh));
    y = OL_g(r_o,4);    % Lambda
    y_bin = zeros(225,2);
    y_bin(:,1) = 0:0.02:4.48;
    for ll = 1:size(y,1)
        for mm = 1:225
            if abs(y(ll) - y_bin(mm)) < 0.02
                y_bin(mm,2) = 1;
              end
        end
    end
    % Clear from empty points
    for pp = 1:225
        if y_bin(pp,2) > 0
            n_bin = n_bin + 1;
            y_bin_edit(n_bin, 1) = Omega_val(hh);
            y_bin_edit(n_bin,2) = y_bin(pp, 1) + 0.01;
        end
    end
end

OL_grain = zeros(size(y_bin_edit, 1), 4);
OL_grain(:,2) = y_bin_edit(:, 1);
OL_grain(:,4) = y_bin_edit(:, 2);

% Make a list cleaned by the projections with too many points
count = histcounts(OL_grain(:,2), 177);
max_num_el = max(count);
n_clean = 0;

% Divide the data in the three regions used for the HT. The regions are designed
% to efficiently cover the different solutions. Division done keeping in mind the % Bragg edges relative to the different hkl values
num_1 = 0;
num_2 = 0;
num_3 = 0;
for aa = 1:size(OL_g, 1)
    if OL_g(aa,4) > 2.5
        num_1 = num_1 + 1;
        OL_g_1(num_1, :) = OL_g(aa, :);
    end
    if OL_g(aa,4) > 2 && OL_g(aa,4) < 3
        num_2 = num_2 + 1;
        OL_g_2(num_2, :) = OL_g(aa, :);
    end
    if OL_g(aa,4) > 1.6 && OL_g(aa,4) < 2.1
        num_3 = num_3 + 1;
        OL_g_3(num_3, :) = OL_g(aa, :);
    end
end

% Find the number of points used in the HT
OL_clean = 15_F_count_points_cleaned(OL_grain, cut_proj);
OL_1 = 15_F_count_points_cleaned(OL_g_1, cut_proj);
OL_2 = 15_F_count_points_cleaned(OL_g_2, cut_proj);
OL_3 = 15_F_count_points_cleaned(OL_g_3, cut_proj);
Number_points = size(OL_clean, 1);

% Calculate the accumulator matrix
%[CM] = Funct_acc_matrix(OL_g, n_large_pixels, steps_alpha, steps_size);
%[CM_clean] = Funct_acc_matrix(OL_clean, n_large_pixels, steps_alpha, steps_size);
[CM1] = 15_Funct_acc_matrix(OL_g_1, n_large_pixels, steps_alpha, steps_size);
[CM2] = 15_Funct_acc_matrix(OL_g_2, n_large_pixels_2, steps_alpha, steps_size);
[CM3] = 15_Funct_acc_matrix(OL_g_3, n_large_pixels_3, steps_alpha, steps_size);

%CM_1 = CM';
%CM_clean_1 = CM_clean';
CM1_1 = CM1';
CM2_1 = CM2';
CM3_1 = CM3';

%g = figure();
%h = pcolor(CM_1); title('Accumulator matrix, all data');
%i = figure();
%l = pcolor(CM_clean_1); title('Accumulator matrix, cleaned data');
m = figure();
n = pcolor(CM1_1); title('Accumulator matrix, l > 2.5');
o = figure();
p = pcolor(CM2_1); title('Accumulator matrix, l in [2,3]');
q = figure();
r = pcolor(CM3_1); title('Accumulator matrix, l in [1.6,2.1]');

% Repeat for l > 2.5 A, l in [1.5, 3] A
% Find the maxima in the accumulator matrix and plot the fitting function
[R_1, alpha_1] = 15_Funct_fit_HT(CM1_1, threshold, n_large_pixels, L_max, ...
    OL_1, threshold_distance_fit);

figure;
plot(OL_g_1(:,2), OL_g_1(:,4), '*'); hold on;
number_centers_1 = size(R_1, 2);
f_1 = zeros(number_centers_1, 180);
for xx = 1:number_centers_1
    for deg = 1:180
        % f is the fitting function
        f_1(xx, deg) = double(R_1(xx)*cosd(alpha_1(xx) + deg));
    end
    plot(1:180, f_1(xx,:), '.'); hold on;
    boundedline(1:180, f_1(xx,:), threshold_distance_fit, 'alpha');
end
xlabel('Omega'); ylabel('Lambda');
axis([0 180 1 4.5]);
title('HT l > 2.5 A, cleaned points');

[R_2, alpha_2] = 15_Funct_fit_HT(CM2_1, threshold_2, n_large_pixels_2, L_max, ...
     OL_2, threshold_distance_fit);

figure;
plot(OL_g_2(:,2), OL_g_2(:,4), '*'); hold on;
number_centers_2 = size(R_2, 2);
f_2 = zeros(number_centers_2, 180);
for xx = 1:number_centers_2
    for deg = 1:180
        % f is the fitting function
        f_2(xx, deg) = double(R_2(xx)*cosd(alpha_2(xx) + deg));
    end
    plot(1:180, f_2(xx,:), '.'); hold on;
    boundedline(1:180, f_2(xx,:), threshold_distance_fit, 'alpha');
end
xlabel('Omega'); ylabel('Lambda');
axis([0 180 1 4.5]);
title('HT l in [2, 3] A, cleaned points');

[R_3, alpha_3] = 15_Funct_fit_HT(CM3_1, threshold_3, n_large_pixels_3, L_max, ...
     OL_3, threshold_distance_fit);

figure;
plot(OL_g_3(:,2), OL_g_3(:,4), '*'); hold on;
number_centers_3 = size(R_3, 2);
f_3 = zeros(number_centers_3, 180);
for xx = 1:number_centers_3
    for deg = 1:180
        % f is the fitting function
        f_3(xx, deg) = double(R_3(xx)*cosd(alpha_3(xx) + deg));
    end
    plot(1:180, f_3(xx,:), '.'); hold on;
    boundedline(1:180, f_3(xx,:), threshold_distance_fit, 'alpha');
end
xlabel('Omega'); ylabel('Lambda');
axis([0 180 1 4.5]);
title('HT l in [1.6, 2.1] A, cleaned points');

% Plot together the results for the two separate lambda ranges
figure;
plot(OL_g(:,2), OL_g(:,4), 'b*'); hold on;
ylim([0 4.5]);
for xx = 1:number_centers_1
    plot(1:180, f_1(xx,:), '.'); hold on;
    boundedline(1:180, f_1(xx,:), threshold_distance_fit, 'alpha');
    hold on;
end
for yy = 1:number_centers_2
    plot(1:180, f_2(yy,:), '.'); hold on;
    boundedline(1:180, f_2(yy,:), threshold_distance_fit, 'alpha');
    hold on;
end
for zz = 1:number_centers_3
    plot(1:180, f_3(zz,:), '.'); hold on;
    boundedline(1:180, f_3(zz,:), threshold_distance_fit, 'alpha');
end
hold on;
plot([0 180 ], [4.054 4.054]);
plot([0 180 ], [2.87 2.87]);
plot([0 180 ], [2.34 2.34]);
plot([0 180 ], [2.03 2.03]);
plot([0 180 ], [1.81 1.81]);
plot([0 180 ], [1.65 1.65]);
xlabel('Omega'); ylabel('Lambda');
axis([0 180 1 4.5]);
title('HT l in [1.6,2.1]; [2,3] and > 2.5 A, cleaned points');
