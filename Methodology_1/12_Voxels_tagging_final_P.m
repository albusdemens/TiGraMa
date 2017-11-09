% Alberto Cereser, 16 October 2015
% Department of Physics, Technical University of Denmark
% alcer@fysik.dtu.dk

% This script returns a voxellized reconstruction of the Fe sample. This is
% done in two steps: at first each voxel is labeled with all the possible
% corresponding grains, and then the most probable one is selected. Version for
% Panda

clear; close all;

center = 127;       % Define the expected rotation centre
min_num_img = 9;    % Requirement on the min number of images for a voxels
dataset_CM = load('3D_reconstr_input.txt');
% Column 1: Projection #; C2: image num; C3: X_cm; C4: Y_cm; C5: tag
num_lines = numel(dataset_CM(:,1));

% Let's make a matrix where to store the voxels
Proj_n = dataset_CM(:,1);
Im_number = dataset_CM(:,2);
X_det = dataset_CM(:,3);
Y_det = dataset_CM(:,4);
Blob_n = dataset_CM(:,5);
x_det_min = min(X_det);
x_det_max = max(X_det);
% diam_sample = max(Y_det(:)) - min(Y_det(:));
n_voxels_x = 200;   % We wanto to keep things in scale, so 1 pixel size =
% size 1 voxel
n_voxels_y = 200;
n_voxels_z = 500;
possible_grains = unique(Blob_n);
n_possible_grains = nnz(possible_grains);
x_center = n_voxels_x/2;
y_center = n_voxels_y/2;

num_el_blob = zeros(n_possible_grains, 1);
% Count the number of angles per grain
for ij = 1:n_possible_grains
    [r, ] = find(Blob_n(:,1) == ij);
    num_el_blob(ij) = nnz(r);
end

Matrix_voxels = zeros(n_voxels_x, n_voxels_y, n_voxels_z);

L = [0. center 0.]'; %[0. center 256/2]; % unit pixels

for i = (1:n_voxels_z)
    disp(i);
    tag = zeros(n_voxels_x, n_voxels_y);
    Matrix_one_z = zeros(n_voxels_x, n_voxels_y, n_possible_grains, 2);
    % Here we store, for each blob, a list of all the possible tag:
    % - a counter, showing how many times the voxel is associated with spots of the
    % considerd grain number;
    % - a percentage, calculated as counter/total_number_images_grain
    % !!! This could be improved by considering what are the possible candidates !!!
    for aa = 1:n_possible_grains
        selected_grain = possible_grains(aa);
        blob_num = 0;
        for j = (1:num_lines)
            dist_points = abs(i-X_det(j));
            if (Blob_n(j) == selected_grain) && dist_points < 100
                Omega = Proj_n(j);
                Spot = Blob_n(j);
                Fit_number = Im_number(j);
                % We check that the distance between the spot CM and the
                % considered voxel is smaller than the extinction spot width
                image = sprintf('/Comb_blobs_final_grouped/Grain_%02i/Comb_final_%03i_%02i.png', aa, Omega, Fit_number);
                IM = imread(image);
                s = size(IM);
                %display(selected_grain);
                for l = 1:n_voxels_x          % unit voxels
                    for k = 1:n_voxels_y      % unit voxels
                        xsam(1) = l/n_voxels_x - 0.5;  % unit cm. The idea is voxels -->
                        % cm --> pixels
                        xsam(2) = k/n_voxels_x - 0.5;  % unit cm
                        xsam(3) = i/n_voxels_x;        % unit cm
                        if (sqrt((xsam(1))^2 + (xsam(2))^2) < 0.5) % The point must be inside the circle
                            O = [cosd(Omega) -sind(Omega) 0.; sind(Omega) cosd(Omega) 0.; 0. 0. 1.];
                            ixdet =  O*xsam'/0.005 + L;  % unit pixels
                            y_det_calc = ixdet(2);
                            z_det_calc = ixdet(3);
                            % We now take into account that the detector image
                            % is rotated by 90 degrees
                            x_det_rot = z_det_calc;
                            y_det_rot = y_det_calc;

                            if (y_det_rot > 0) && (y_det_rot < s(1))
                                if IM(floor(y_det_rot), floor(x_det_rot)) > 0
                                    %blob_num = selected_grain;
                                    Matrix_one_z(k, l, aa, 1) = Matrix_one_z(k, l, aa, 1) + 1;
                                end
                            end
                        end
                    end
                end
            end
        end
    end
    % We want to limit our analysis to when the voxels have a label
    for ll = 1:n_voxels_x          % unit voxels
        for kk = 1:n_voxels_y      % unit voxels
            if nnz(Matrix_one_z(ll, kk, :, 1)) > 0 % We need at least one possible tag
                for jj = 1:n_possible_grains
                    if Matrix_one_z(ll, kk, jj, 1) > num_el_blob(jj)*0.6
                        Matrix_one_z(ll, kk, jj, 2) = Matrix_one_z(ll, kk, jj, 1);
                    else
                        Matrix_one_z(ll, kk, jj, 2) = 0;
                    end
                end
                if nnz(Matrix_one_z(ll, kk, :, 2)) > 0
                    [possible_tag, ind_possible_tag] = max(Matrix_one_z(ll, kk, :, 2));
                    tag(ll,kk) = ind_possible_tag;
                end
            end
        end
    end
    IM_layer = sprintf('/Grain_reconstruction/Grain_reconstruction_NN/Final_shape_Fe_%03i.fits', i);
    fitswrite(tag, IM_layer);
end
