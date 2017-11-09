% Alberto Cereser, 30 November 2015
% Technical University of Denmark, alcer@fysik.dtu.dk

% Script to count the number of grains in the reconstruction obtained
% combining the [100;500], [100;1000], [100;2000], 1000+ datasets.
% Approach: load the biggest blobs, find which of the smaller one are inside the
% bigger ones, iterate.

close all; clear;

grain_num = 17; % Possible to save vtk file relative to a single grain

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

%savevtk(V1_refined,'Numbered_grains.vtk');
disp('V1 done');

unq = unique(V1_refined);
for nn = 2:size(unique(V1_refined)) % First value is 0
    [xx1, yy1, zz1] = ind2sub(size(V1_refined),find(V1_refined == unq(nn)));
    x_cm1(nn-1) = mean(xx1);
    y_cm1(nn-1) = mean(yy1);
    z_cm1(nn-1) = mean(zz1);
end

figure; scatter3(x_cm1', y_cm1', z_cm1');
CM_grains = [x_cm1' y_cm1' z_cm1'] ;

%%%% Write Cm grains
%dlmwrite('CM_grains.txt', CM_grains, 'delimiter', '\t', 'precision', 6);

% Calculate the volume of each grain
List_labels = unique(V1_refined);
num_labels = nnz(List_labels);
Vol_counter = zeros(num_labels, 2);
for uu = 1:num_labels
    lab_num = List_labels(uu+1);
    num_vol = 0;
    for vv = 1:200
        for ww = 1:200
            for zz = 1:n_1000_plus-3
                if V1_refined(vv, ww, zz) == lab_num
                    num_vol = num_vol + 1;
                end
            end
        end
    end
    Vol_counter(uu, 1) = lab_num;
    Vol_counter(uu, 2) = num_vol;
end

%%%% Write volume grains
%dlmwrite('Volume_grains.txt', Vol_counter, 'delimiter', '\t');

% Mark grains using progressive numbering
Gr_tag = unique(V1);
V2 = zeros(size(V1));
for hh = 2:size(Gr_tag) % First value is 0
    Gr_n = Gr_tag(hh);
    for ii = 1:size(V1,1)
        for jj = 1:size(V1,2)
            for kk = 1:size(V1,3)
                if V1(ii,jj,kk) == Gr_n
                    V2(ii,jj,kk) = hh;
                end
            end
        end
    end
end
%savevtk(V2,'Grains_rescaled.vtk');

% Select a certain grain
[X,Y,Z] = size(V1);
Selected_G = zeros(size(V1));
for aa = 1:X
    for bb = 1:Y
        for cc = 1:Z
            if V2(aa, bb, cc) == grain_num
             Selected_G(aa, bb, cc) = 2;%V2(aa, bb, cc);
            end
        end
    end
end

%%%% Save reconstruction for a certain grain
%savevtk(Selected_G,'G_old_reconstruction.vtk');
[Info_CM] = Function_CM_grains(Selected_G)
return;

% Search what is the final ID for each grain of the various I intervals
grains_final = unique(V1);
n_gfinal = nnz(grains_final);
List_grain_fin = zeros(n_gfinal, 100);
List_grain_fin(:,1) = grains_final(2:size(grains_final, 1), 1);
Grain_size = zeros(n_gfinal,2);
Grain_size(:,1) = grains_final(2:size(grains_final,1));
for ii = 2:n_gfinal
    g_elements = 0;   % To calculate the grain volume, count the contained voxels
    for aa = 1:size(Vol_100_500, 2)
        for bb = 1:size(Vol_100_500, 2)
            for cc = 1:n_1000_plus-3
                if V1(aa, bb, cc) == grains_final(ii)
                    g_elements = g_elements + 1;
                    if Vol_1000_plus(aa, bb, cc) > 0
                        T_1 = Vol_1000_plus(aa, bb, cc);
                    end
                    if Vol_100_2000(aa, bb, cc) > 0
                        T_2 = Vol_100_2000(aa, bb, cc);
                    end
                    if Vol_100_1000(aa, bb, cc) > 0
                        T_3 = Vol_100_1000(aa, bb, cc);
                    end
                    if Vol_100_500(aa, bb, cc) > 0
                        T_4 = Vol_100_500(aa, bb, cc);
                    end
                    tags_grain = List_grain_fin(ii-1, 2:size(List_grain_fin, 2));
                    if Vol_1000_plus(aa, bb, cc) > 0 && isempty(find(tags_grain == T_1)) == 1
                        tags_already = nnz(List_grain_fin(ii-1,:));
                        List_grain_fin(ii-1, tags_already + 1) = T_1;
                    end
                    if Vol_100_2000(aa, bb, cc) > 0 && isempty(find(tags_grain == T_2)) == 1
                        tags_already = nnz(List_grain_fin(ii-1,:));
                        List_grain_fin(ii-1, tags_already + 1) = T_2;
                    end
                    if Vol_100_1000(aa, bb, cc) > 0 && isempty(find(tags_grain == T_3)) == 1
                        tags_already = nnz(List_grain_fin(ii-1,:));
                        List_grain_fin(ii-1, tags_already + 1) = T_3;
                    end
                    if Vol_100_500(aa, bb, cc) > 0 && isempty(find(tags_grain == T_4)) == 1
                        tags_already = nnz(List_grain_fin(ii-1,:));
                        List_grain_fin(ii-1, tags_already + 1) = T_4;
                    end
                end
            end
        end
    end
    Grain_size(ii, 2) = g_elements;
end
% Clean List_grain_fin and save matrix
List_grain_fin(:, ~any(List_grain_fin,1) ) = [];
% First column: final ID of the grain; following cols: ID in the partial
% reconstructions
dlmwrite('Grains_size.txt', Grain_size, 'delimiter', '\t');
dlmwrite('Grains_combination.txt', List_grain_fin, 'delimiter', '\t');

% Unify OL lists using cat
A_1 = load('Reconstruction_1000/Omega_lambda_final.txt');
A_2 = load('Reconstruction_100/Omega_lambda_final.txt');
A_3 = load('Reconstruction_100-1000/Omega_lambda_final.txt');
A_4 = load('Reconstruction_100-500/Omega_lambda_final.txt');

delete('OL_final_grains/*');
for ii = 1:size(List_grain_fin(:,1))
    OL_grain = zeros(1, 6);
    t_grain = List_grain_fin(ii, 2:size(List_grain_fin,2));
    num_t = nnz(t_grain);
    for jj = 1:num_t
        if t_grain(jj) <= m_1
            [r1, ] = find(A_1(:,1) == t_grain(jj));
            Data_grain1 = zeros(size(r1,1), 6);
            Data_grain1(1:size(r1,1),1:4) = A_1(r1, 1:4);
            Data_grain1(1:size(r1,1),5) = List_grain_fin(ii,1);
            Data_grain1(1:size(r1,1),6) = 1;
            OL_grain = cat(1, OL_grain, Data_grain1);
        end
        if t_grain(jj) > m_1 && t_grain(jj) <= (m_1 + m_2)
            [r2, ] = find(A_2(:,1) == (t_grain(jj) - m_1));
            Data_grain2 = zeros(size(r2,1), 6);
            Data_grain2(1:size(r2,1),1:4) = A_2(r2, 1:4);
            Data_grain2(1:size(r2,1),5) = List_grain_fin(ii,1);
            Data_grain2(1:size(r2,1),6) = 2;
            OL_grain = cat(1, OL_grain, Data_grain2);
        end
        if t_grain(jj) > (m_1 + m_2) && t_grain(jj) <= (m_1 + m_2 + m_3)
            [r3, ] = find(A_3(:,1) == (t_grain(jj) - (m_1 + m_2)));
            Data_grain3 = zeros(size(r3,1), 6);
            Data_grain3(1:size(r3,1),1:4) = A_3(r3, 1:4);
            Data_grain3(1:size(r3,1),5) = List_grain_fin(ii,1);
            Data_grain3(1:size(r3,1),6) = 3;
            OL_grain = cat(1, OL_grain, Data_grain3);
        end
        if t_grain(jj) > (m_1 + m_2 + m_3)
            [r4, ] = find(A_4(:,1) == (t_grain(jj) - (m_1 + m_2 + m_3)));
            Data_grain4 = zeros(size(r4,1), 6);
            Data_grain4(1:size(r4,1),1:4) = A_4(r4, 1:4);
            Data_grain4(1:size(r4,1),5) = List_grain_fin(ii,1);
            Data_grain4(1:size(r4,1),6) = 4;
            OL_grain = cat(1, OL_grain, Data_grain4);
        end
    end
    Output_g = sprintf('OL_final_grains/OL_%02i.txt', ii);
    dlmwrite(Output_g, OL_grain, 'delimiter', '\t');
    % C1: Grain number in partial reconstruction; C2: Omega; C3: image
    % number; C4: lambda; C5: grain ID in final reconstruction; C6:
    % partial reconstruction number
end
return;

% Slice the sample at different angles. To do so, transform all the voxel
% coordinates in polar values
Polar_voxel = zeros(499*200*200, 6);
point_num = 0;
for kk = 1:499
    for ii = 1:200
        for jj = 1:200
            point_num = point_num + 1;
            Polar_voxel(point_num, 1:3) = [ii jj kk];
            [alpha, rho] = cart2pol(ii - 100, jj - 100);
            Polar_voxel(point_num, 4:5) = [radtodeg(alpha) rho];
            Polar_voxel(point_num, 6) = V1(ii,jj,kk);
        end
    end
end

Coord_voxel = zeros(499*200*200, 4);
voxel_n = 0;
for ii = 1:200
    for jj = 1:200
        for kk = 1:499
            if V1(ii,jj,kk) > 0
                voxel_n = voxel_n + 1;
                Coord_voxel(voxel_n, 1:3) = [ii, jj, kk];
                Coord_voxel(voxel_n, 4) = V1(ii,jj,kk);
            end
        end
    end
end

%for ii = 1:200
%    for jj = 1:200
%        for kk = size(V1,3)
%            if V1(ii,jj,kk) == 0
%                V2(ii,jj,kk) = NaN;
%            else
%                V2(ii,jj,kk) = V1(ii,jj,kk);
%            end
%        end
%    end
%end

% figure;
% cmap = colormap; % get current colormap
% cmap=cmap([1 2],:); % set your range here
% colormap(cmap); % apply new colormap
% hpat = PATCH_3Darray(V1, 'col');
% colorbar;

return;

% Save the slices
for ij = 1:1:180
    Slice = zeros(200, 500);
    for jk = 1:size(Polar_voxel, 1)
        if round(Polar_voxel(jk,4)) == ij && Polar_voxel(jk,5) < 99
            if Polar_voxel(jk, 6) > 0;
                Slice(round(Polar_voxel(jk,5)) + 100, Polar_voxel(jk,3)) = Polar_voxel(jk, 6);
            end
        elseif round(Polar_voxel(jk,4)) == ij - 180 && abs(Polar_voxel(jk,5)) < 99
            if Polar_voxel(jk, 6) > 0;
                Slice(-round(Polar_voxel(jk,5)) + 100, Polar_voxel(jk,3)) = Polar_voxel(jk, 6);
            end
        end
    end
    filename = sprintf('Slices/Slice_%03i.fits', ij);
    fitswrite(Slice', filename);
end
