% Alberto Cereser, 2 October 2015
% Department of Physics, Technical University of Denmark
% alcer@fysik.dtu.dk

% Script to plot the CM of the extinction spots and to group them. The
% grouping procedure is designed to combine CM values from the Hough
% transform that have similar X and Y component, but different Z

%clear; close all;

function 11_Plot_and_group_CM_Murofi(a_ell, b_ell, c_ell)

min_n_points = 3;  % Min number of points to be considered to calculate
% the CM of the voxel cluster

% list_CM_Murofi = sprintf('/Users/Alberto/Documents/Data_analysis/J-PARC_2015/Data_analysis_final/CM_grouped/List_CM_all.txt');
% List_spots = load(list_CM_Murofi);
%
% X = List_spots(:,1);
% Y = List_spots(:,2);
% Omega = List_spots(:,4);
%
% figure; scatter3(X, Y, Omega); title 'CM of the spots'; xlabel 'Y'; ylabel 'X'; zlabel 'Omega';

List_spots_final = load('CM_alpha_R.txt');

%Remove the zero rows
List_spots_final( ~any(List_spots_final,2), : ) = [];

Z = List_spots_final(:,4);
Angle = List_spots_final(:,6);
R = List_spots_final(:,7);

for i = 1:size(Z)
    X_final(i) = R(i)*(cosd(Angle(i)));
    Y_final(i) = R(i)*(sind(Angle(i)));
end

% For all numbers, the scale is mm
figure; scatter3(X_final, Y_final, Z); title 'CM of the spots'; xlabel 'Y'; ylabel 'X'; zlabel 'Z';
%figure; scatter(X_final, Y_final);

% Plot a 2D accumulator matrix
X_acc = (X_final + 2.5)*20;
Y_acc = (Y_final + 2.55)*20;
Coord_final = [X_acc', Y_acc'];

Count_matrix = zeros(100,100);
for hh = 1:size(Coord_final, 1)
    x_point = round(Coord_final(hh,1));
    y_point = round(Coord_final(hh,2));
    Count_matrix(x_point, y_point) = Count_matrix(x_point, y_point) + 1;
end

% Plot the accumulation matrix
figure; mesh(Count_matrix);

% 3D accumulation matrix
Count_matrix_3D = zeros(100,100,100);
for hh = 1:size(Coord_final, 1)
    x_point = round(Coord_final(hh,1));
    y_point = round(Coord_final(hh,2));
    z_point = round(Z(hh));
    Count_matrix_3D(x_point, y_point, z_point) = 1;
end

% Use an ellipsoidal object to find the CM of the point clusters
[X_one, Y_one, Z_one] = ind2sub(size(Count_matrix_3D),find(Count_matrix_3D == 1));
num_ones = nnz(X_one);
% Counter is the matrix where we store, for each nnz value, the number of
% neighbours in the considered ellipsoid
Counter = zeros(num_ones, 4);
Counter(:,1:3) = [X_one, Y_one, Z_one];
for aa = 1:num_ones
    for bb = 1:num_ones
        if aa ~= bb
           dist_2 = ((X_one(aa) - X_one(bb))/a_ell)^2 + ((Y_one(aa) - ...
               Y_one(bb))/b_ell)^2 + ((Z_one(aa) - Z_one(bb))/c_ell)^2;
           if dist_2 < 1
                Counter(aa, 4) = Counter(aa, 4) + 1;
           end
        end
    end
end

% Make a color index for connection level
c = Counter(:,4);

[x_acc_3D, y_acc_3D, z_acc_3D] = ind2sub(size(Count_matrix_3D), find(Count_matrix_3D > 0));
% Plot the calculated CM values together with the raw data
%figure;
%scatter3(x_acc_3D, y_acc_3D, z_acc_3D, 30, c, '*'); hold on;
%colorbar; title('Cm values from the HT and their connectivity');
%xlabel('X'); ylabel('Y'); zlabel('Z'); hold on;
%scatter3(CM(:,1), CM(:,2), CM(:,3), 'o');
% % Plot, for all voxels, the number of elements per ellipsoid centered on them
%figure; histogram(Counter(:, 4));
%title('Number of elements connected to a voxel');
%xlabel('Voxel connection');
%ylabel('Counter');

% Locate the voxels that are at least N-connected (paremeter defined at
% the beginning of the script)
idx_8_conn = find(Counter(:,4) >= min_n_points);
Counter_conn = Counter(idx_8_conn, :);
x_conn = Counter_conn(:,1);
y_conn = Counter_conn(:,2);
z_conn = Counter_conn(:,3);
c_conn = Counter_conn(:,4);

% Plot the 3D distribution of the voxels with at least the selected number
% of connections
figure; scatter3(x_conn, y_conn, z_conn, 30, c_conn, '*'); colorbar; hold on;
title('CM values with at leats N connections');
xlabel('X'); ylabel('Y'); zlabel('Z');

Tag_matrix = zeros(size(Counter_conn,1), 5);
Tag_matrix(:,1:4) = Counter_conn(:,1:4);
% Find the center of mass of each cluster of voxels
for ij = 1:size(idx_8_conn,1)
    for jk = 1:size(idx_8_conn,1)
        if ij ~= jk
           dist_2 = ((x_conn(ij) - x_conn(jk))/a_ell)^2 + ((y_conn(ij) - ...
               y_conn(jk))/b_ell)^2 + ((z_conn(ij) - z_conn(jk))/c_ell)^2;
           if dist_2 < 1
               if Tag_matrix(ij,5) == 0 && Tag_matrix(jk,5) == 0
                   Tag_matrix(ij,5) = ij;
                   Tag_matrix(jk,5) = ij;
               elseif Tag_matrix(ij,5) ~= 0 && Tag_matrix(jk,5) == 0
                   Tag_matrix(jk,5) = Tag_matrix(ij,5);
               elseif Tag_matrix(ij,5) == 0 && Tag_matrix(jk,5) ~= 0
                   Tag_matrix(ij,5) = Tag_matrix(jk,5);
               end
           end
        end
    end
end

Cluster_tags = unique(Tag_matrix(:,5));
Cluster_tags(~any(Cluster_tags,2), : ) = [];    % Delete zero row
Num_tagged_clusters = nnz(Cluster_tags);
% For each cluster, find the CM of the voxels with at least 15 neighbours
Clusters_CM = zeros(Num_tagged_clusters, 4);
for i = 1:Num_tagged_clusters
    [r,~] = find(Tag_matrix(:,5) == Cluster_tags(i));
    X_CM = mean(x_conn(r));
    Y_CM = mean(y_conn(r));
    Z_CM = mean(z_conn(r));
    Clusters_CM(i, :) = [Cluster_tags(i), X_CM, Y_CM, Z_CM];
end

% For each CM, we build an ellipsoid and calculate a new CM, this time
% taking into account all the initial data
All_data_tagged = zeros(num_ones, 4);
for i = 1:Num_tagged_clusters
    for ii = 1:num_ones
        dist_3 = ((Clusters_CM(i,2) - x_acc_3D(ii))/a_ell)^2 + ((Clusters_CM(i,3) - ...
            y_acc_3D(ii))/b_ell)^2 + ((Clusters_CM(i,4) - z_acc_3D(ii))/c_ell)^2;
        All_data_tagged(ii, 1:3) = Counter(ii, 2:4);
        if dist_3 < 1
            All_data_tagged(ii, 4) = i;
        end
    end
end

% Plot all the data, with the points composing the clusters in colorscale
c_tags_n = All_data_tagged(:,4);
% figure; scatter3(x_acc_3D, y_acc_3D, z_acc_3D, 30, c_tags_n, '*'); colorbar; hold on;
% title('All points, with clusters in colorscale');
% xlabel('X'); ylabel('Y'); zlabel('Z');

% For each cluster, calculate the CM using all the available points
for i = 1:Num_tagged_clusters
    [r,~] = find(All_data_tagged(:,4) == i);
    X_CM_fin = mean(x_acc_3D(r));
    Y_CM_fin = mean(y_acc_3D(r));
    Z_CM_fin = mean(z_acc_3D(r));
    Clusters_CM_fin(i, :) = [i, X_CM_fin, Y_CM_fin, Z_CM_fin];
end

% Plot all data together with the definitive CM values
% Visually check the tagged clusters distribution
c_tags = Tag_matrix(:,5);
figure; scatter3(x_conn, y_conn, z_conn, 30, c_tags, '*'); colorbar; hold on;
% title('Distribution of the tagged CM values');
xlabel('X'); ylabel('Y'); zlabel('Z'); hold on;
scatter3(Clusters_CM_fin(:,2), Clusters_CM_fin(:,3), Clusters_CM_fin(:,4), 'o');

% Plot the CM alone
figure; scatter3(Clusters_CM_fin(:,2), Clusters_CM_fin(:,3), Clusters_CM_fin(:,4), 'o');
title('Distribution of the CM values'); xlabel('X'); ylabel('Y'); zlabel('Z');

% Write CM of the clusters on a file
dlmwrite('CM_clusters_final.txt', Clusters_CM_fin, 'delimiter',' ', 'precision', 6);
