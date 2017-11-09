% Alberto Cereser, 13 October 2015
% Technical University of Denmark, alcer@fysik.dtu.dk

% Script that cleans the output of the HT, made by Hough_trans_pol_rolling_final.m
% Aims, by listing the images relative to a solution. For certain
% projections, there are more than one image. These are combined by
% Comb_im_before_reconstr.m

function 11_Clean_result_HT(a_ell, b_ell, c_ell)

Data_HT = load('CM_alpha_R.txt');
X_val = unique(Data_HT(:,2));
X_val_edit = zeros(size(X_val, 1) - 1, 1);
line_edit = 0;
for i = 1:size(X_val,1)
    if X_val(i,1) > 0
        line_edit = line_edit + 1;
        X_val_edit(line_edit,1) = X_val(i,1);
    end
end

% For each rolling interval, find the blobs relative to each solution
point_num = 0;
sol_num = 0;
for ii = 1:size(X_val_edit,1)
    % Find data relative to the considered interval
    [r_int, ] = find(Data_HT(:, 2) == X_val_edit(ii,1));
    % Find the unique (alpha, R combinations)
    Comb = unique(Data_HT(r_int, 6:7), 'rows');
    % Store data relative to the considered interval in a separate matrix
    Data_HT_roll = zeros(size(r_int, 1),7);
    Data_HT_roll(:,:) = Data_HT(r_int, :);
    for jj = 1:size(Comb, 1)
        sol_num = sol_num + 1;
        [r, ] = find(Data_HT_roll(:, 6) == Comb(jj, 1) & Data_HT_roll(:, 7) == Comb(jj, 2));
        for kk = 1:size(r,1)
           point_num = point_num + 1;
           % Stored data: interval center, line num in Data_HT, Alpha and R
           % from the HT, omega and image number for the contributor
           HT_Elements(point_num, 1:6) = [X_val_edit(ii,1), sol_num, ...
               Comb(jj, 1:2), Data_HT_roll(r(kk), 1), Data_HT_roll(r(kk), 5)];
        end
    end
end

% Load the values of the CM, grouped by Plot_and_group_CM_Murofi.m, and
% return the data relative to the components
CM_final = load('CM_clusters_final.txt');
Z_scaled = HT_Elements(:,1);
Angle = HT_Elements(:,3);
R = HT_Elements(:,4);
X_final = zeros(size(HT_Elements, 1), 1);
Y_final = zeros(size(HT_Elements, 1), 1);
for ab = 1:size(HT_Elements, 1)
    X_final(ab) = R(ab)*(cosd(Angle(ab)));
    Y_final(ab) = R(ab)*(sind(Angle(ab)));
end
X_scaled = (X_final + 2.5)*20;
Y_scaled = (Y_final + 2.5)*20;

num_el_fin = 0;
for aa = 1:size(CM_final,1)
    for bb = 1:size(HT_Elements, 1)
        dist = ((CM_final(aa,2) - X_scaled(bb))/a_ell)^2 + ((CM_final(aa,3) - ...
            Y_scaled(bb))/b_ell)^2 + ((CM_final(aa,4) - Z_scaled(bb))/c_ell)^2;
        if dist < 1
            num_el_fin = num_el_fin + 1;
            HT_Elements_fin(num_el_fin, 1:6) = HT_Elements(bb, 1:6);
            % Columns: X coord, solution number (from HT), alpha, R, omega,
            % image ID, grain number
            HT_Elements_fin(num_el_fin, 7) = aa;    % Grain number
        end
    end
end

% Plot data in HT_Elements_fin to check that everything is fine
Z_fin = HT_Elements_fin(:,1);
Angle_fin = HT_Elements_fin(:,3);
R_fin = HT_Elements_fin(:,4);
X_final = zeros(size(HT_Elements_fin, 1), 1);
Y_final_fin = zeros(size(HT_Elements_fin, 1), 1);
for ab = 1:size(HT_Elements_fin, 1)
    X_final_fin(ab) = R_fin(ab)*(cosd(Angle_fin(ab)));
    Y_final_fin(ab) = R_fin(ab)*(sind(Angle_fin(ab)));
end
X_scaled_fin = (X_final_fin + 2.5)*20;
Y_scaled_fin = (Y_final_fin + 2.5)*20;
%figure; scatter3(X_scaled_fin, Y_scaled_fin, Z_fin);

% Save HT_Elements
dlmwrite('Cleaned_results_HT.txt', HT_Elements_fin, 'delimiter',' ', 'precision', 7);
