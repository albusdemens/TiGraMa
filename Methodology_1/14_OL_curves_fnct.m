% Alberto Cereser, 11 May 2016
% Technical University of Denmark, alcer@fysik.dtu.dk

% Script to count the number of grains in the reconstruction obtained
% combining the [100;500], [100;1000], [100;2000], 1000+ datasets.
% Approach: load the biggest blobs, find which of the smaller one are inside the
% bigger ones, iterate.

% Based on Count_grains.m, called by Find_grains_ID.m

function [V1_refined, V, Vol_1000_plus, Vol_100_2000, Vol_100_1000, Vol_100_500] = 14_OL_curves_fnct

% I > 1000
img_1000_plus = dir(sprintf('Reconstruction_1000/Grain_reconstruction/Reconstruction_0_59/'));
n_1000_plus = length(img_1000_plus);    % Number of layers
path_1000_plus = sprintf('Reconstruction_1000/Grain_reconstruction/Reconstruction_0_59/');

for aa = 1:n_1000_plus-3
    currentfilename = img_1000_plus(aa+3).name;
    path_image = strcat(path_1000_plus, currentfilename);
    A = double(fitsread(path_image));
    Vol_1000_plus(:,:,aa) = A(:,:);
end
clear A;

% I in [100, 2000]
img_100_2000 = dir(sprintf('Reconstruction_100/Grains_reconstruction/'));
path_100_2000 = sprintf('Reconstruction_100/Grains_reconstruction/');

add_1 = size(unique(nnz(Vol_1000_plus)),1) + 1;
for aa = 1:n_1000_plus-3
    currentfilename = img_100_2000(aa+3).name;
    path_image = strcat(path_100_2000, currentfilename);
    A = double(fitsread(path_image));
    for uu = 1:size(A,1)
        for vv = 1:size(A,2)
            if A (uu, vv) > 0
                Vol_100_2000(uu,vv,aa) = A(uu,vv) + add_1;
            else
                Vol_100_2000(uu,vv,aa) = 0;
            end
        end
    end
end
clear A;

% I in [100, 1000]
img_100_1000 = dir(sprintf('Reconstruction_100-1000/Grains_reconstruction'));
path_100_1000 = sprintf('Reconstruction_100-1000/Grains_reconstruction/');

add_2 = size(unique(nnz(Vol_1000_plus)),1) + size(unique(nnz(Vol_100_2000)),1) + 2;
for aa = 1:n_1000_plus-3
    currentfilename = img_100_1000(aa+3).name;
    path_image = strcat(path_100_1000, currentfilename);
    A = double(fitsread(path_image));
    for uu = 1:size(A,1)
        for vv = 1:size(A,2)
            if A (uu, vv) > 0
                Vol_100_1000(uu,vv,aa) = A(uu,vv) + add_2;
            else
                Vol_100_1000(uu,vv,aa) = 0;
            end
        end
    end
end
clear A;

% I in [100, 500]
img_100_500 = dir(sprintf('Reconstruction_100-500/Grains_reconstruction/'));
path_100_500 = sprintf('Reconstruction_100-500/Grains_reconstruction/');

add_3 = size(unique(nnz(Vol_1000_plus)),1) + size(nnz(unique(Vol_100_2000)),1) + size(unique(nnz(Vol_100_1000)),1) + 3;
for aa = 1:n_1000_plus-3
    currentfilename = img_100_500(aa+3).name;
    path_image = strcat(path_100_500, currentfilename);
    A = double(fitsread(path_image));
    for uu = 1:size(A,1)
        for vv = 1:size(A,2)
            if A (uu, vv) > 0
                Vol_100_500(uu,vv,aa) = A(uu,vv) + add_3;
            else
                Vol_100_500(uu,vv,aa) = 0;
            end
        end
    end
end
clear A;

% Incrementally, find grains of the smaller datasets that are outisde the bigger
% ones
V_comb_1 = zeros(size(Vol_100_500));
V_comb_2 = zeros(size(Vol_100_500));
V_comb_3 = zeros(size(Vol_100_500));
V = zeros(size(Vol_100_500));
V1 = zeros(size(Vol_100_500));
V2 = zeros(size(Vol_100_500));

Max_V0 = max(max(max(Vol_1000_plus))) + 1;
for ii = 1:size(Vol_100_500, 2)
    for jj = 1:size(Vol_100_500, 2)
        for kk = 1:n_1000_plus-3
            if (Vol_1000_plus(ii, jj, kk) > 0 && Vol_100_2000(ii,jj,kk) > 0)
                V_comb_1(ii,jj,kk) = Vol_1000_plus(ii, jj, kk);
            elseif (Vol_1000_plus(ii, jj, kk) > 0 && Vol_100_2000(ii,jj,kk) == 0)
                V_comb_1(ii,jj,kk) = Vol_1000_plus(ii, jj, kk);
            elseif (Vol_1000_plus(ii, jj, kk) == 0 && Vol_100_2000(ii,jj,kk) > 0)
                V_comb_1(ii,jj,kk) = Vol_100_2000(ii, jj, kk) + Max_V0;
            end
        end
    end
end

Max_V1 = max(max(max(V_comb_1))) + 1;
for ii = 1:size(Vol_100_500, 2)
    for jj = 1:size(Vol_100_500, 2)
        for kk = 1:n_1000_plus-3
            if (V_comb_1(ii, jj, kk) > 0 && Vol_100_1000(ii,jj,kk) > 0)
                V_comb_2(ii,jj,kk) = V_comb_1(ii, jj, kk);
            elseif (V_comb_1(ii, jj, kk) > 0 && Vol_100_1000(ii,jj,kk) == 0)
                V_comb_2(ii,jj,kk) = V_comb_1(ii, jj, kk);
            elseif (V_comb_1(ii, jj, kk) == 0 && Vol_100_1000(ii,jj,kk) > 0)
                V_comb_2(ii,jj,kk) = Vol_100_1000(ii, jj, kk) + Max_V1;
            end
        end
    end
end
Max_V2 = max(max(max(V_comb_2))) + 1;

for ii = 1:size(Vol_100_500, 2)
    for jj = 1:size(Vol_100_500, 2)
        for kk = 1:n_1000_plus-3
            if (V_comb_2(ii, jj, kk) > 0 && Vol_100_500(ii,jj,kk) > 0)
                V_comb_3(ii,jj,kk) = V_comb_2(ii, jj, kk);
            elseif (V_comb_2(ii, jj, kk) > 0 && Vol_100_500(ii,jj,kk) == 0)
                V_comb_3(ii,jj,kk) = V_comb_2(ii, jj, kk);
            elseif (V_comb_2(ii, jj, kk) == 0 && Vol_100_500(ii,jj,kk) > 0)
                V_comb_3(ii,jj,kk) = Vol_100_500(ii, jj, kk) + Max_V2;
            end
        end
    end
end

% Introduce V
for ii = 1:200
    for jj = 1:200
        for kk = 1:n_1000_plus-3
            if V_comb_3(ii,jj,kk) == 0
                V(ii,jj,kk) = 0; % NAN;
            else
                V(ii,jj,kk) = V_comb_3(ii,jj,kk);
            end
        end
    end
end

disp('V done');

grains_id = unique(V);
num_grains = size(grains_id,1);
x_cm = zeros(num_grains,1);
y_cm = zeros(num_grains,1);
z_cm = zeros(num_grains,1);

% Assign to a grain the value of its CM
for ll = 1:num_grains
    idx = grains_id(ll);
    xx = 0;
    yy = 0;
    zz = 0;
    [xx, yy, zz] = ind2sub(size(V),find(V == idx));
    x_cm = mean(xx);
    y_cm = mean(yy);
    z_cm = mean(zz);
    cm_ID = V(round(x_cm), round(y_cm), round(z_cm));
    num_vox = size(find(V == idx),1);
    if num_vox > 110
        for dd = 1:size(V,1)
            for ee = 1:size(V,2)
                for ff = 1:size(V,3)
                    if V(dd, ee, ff) == idx
                        V1(dd, ee, ff) = cm_ID;
                    end
                end
            end
        end
    end
end

V1_refined = zeros(size(V1));
G_id = unique(V1);
for cc = 2:size(G_id,1);
    for dd = 1:size(V1,1)
        for ee = 1:size(V1,2)
            for ff = 1:size(V1,3)
                if V1(dd, ee, ff) == G_id(cc)
                    V1_refined(dd, ee, ff) = cc - 1;
                end
            end
        end
    end
end
