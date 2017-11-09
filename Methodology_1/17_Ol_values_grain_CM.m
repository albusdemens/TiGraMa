% Alberto Cereser, 3 January 2016
% Technical University of Denmark, alcer@fysik.dtu.dk

% Script to clean the OL values for a grain using its CM, returned by
% Combine_all_3D_sets.m. Approach: a) for each OL combination, load the
% corresponding image; b) check if the image contains the CM; c) keep only the
% OL values whose corresponding image contains the CM.

close all; clear;
IN_grains = load('CM_grains.txt');
CM_grains(:,1:3) = IN_grains(:,:);
CM_grains(:,4) = 1:size(IN_grains, 1);
OL = load('OL_comb_after_CM.txt');
Grains_ID = CM_grains(:,4);

% Load files listing all the isolated blobs, for each I range
Info_1000_plus = load('Properties_isolated_blobs_1000.txt');
Info_100_2000 = load('Properties_isolated_blobs_100.txt');
Info_100_1000 = load('Properties_isolated_blobs_100-1000.txt');
Info_100_500 = load('Properties_isolated_blobs_100-500.txt');

num_ok = 0;
for ii = 1:size(Grains_ID,1)
    disp(ii/size(Grains_ID,1));
    grain_num = Grains_ID(ii);
    x_grain = CM_grains(ii, 1);
    y_grain = CM_grains(ii, 2);
    z_grain = CM_grains(ii, 3);
    R = sqrt((x_grain-100)^2 + (y_grain-100)^2);
    alpha = acosd((x_grain - 100)/R);
    for jj = 1:size(OL,1)
        if OL(jj, 1) == grain_num
            Omega = OL(jj, 3);
            Image_n = OL(jj, 4);
            Grain_range = OL(jj, 2);    % Depends on the extinction spots area
            if Grain_range == 4
                [r, ] = find((Info_1000_plus(:,1)*3 == Omega) & (Info_1000_plus(:,2) == Image_n));
            elseif Grain_range == 3
                [r, ] = find((Info_100_2000(:,1)*3 == Omega) & (Info_100_2000(:,2) == Image_n));
            elseif Grain_range == 2
                [r, ] = find((Info_100_1000(:,1)*3 == Omega) & (Info_100_1000(:,2) == Image_n));
            elseif Grain_range == 1
                [r, ] = find((Info_100_500(:,1)*3 == Omega) & (Info_100_500(:,2) == Image_n));
            end
            for kk = 1:size(r,1)
                proj = Omega/3;
                if Grain_range == 4
                    file_name = sprintf('Isolated_blobs_1000_all/Isolated_blob_%03i_%05i_%02i.png', proj, Image_n, kk);
                elseif Grain_range == 3
                    file_name = sprintf('Isolated_blobs_100_all/Isolated_blob_%03i_%05i_%02i.png', proj, Image_n, kk);
                elseif Grain_range == 2
                  file_name = sprintf('Isolated_blobs_100-1000_all/Isolated_blob_%03i_%05i_%02i.png', proj, Image_n, kk);
                elseif Grain_range == 1
                    file_name = sprintf('Isolated_blobs_100-500_all/Isolated_blob_%03i_%05i_%02i.png', proj, Image_n, kk);
                end
                IM = imread(file_name);
                x_img = z_grain;
                y_img = R*cosd(alpha+Omega) + 100;
                if IM(round(y_img), round(x_img)) > 0
                    num_ok = num_ok + 1
                    OL_OK(num_ok, :) = OL(jj, :);
                end
            end
        end
    end
end

dlmwrite('OL_final_OK.txt', OL_OK, 'delimiter', '\t');
