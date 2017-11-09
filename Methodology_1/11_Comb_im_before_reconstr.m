% Alberto Cereser, 14 October 2015
% Technical University of Denmark, alcer@fysik.dtu.dk

% In Cleaned_results_HT.txt, for some projections we have more than one
% image. This would result in anisotropies in the reconstruction (certain
% are more important than others). To solve this, we combine spots before
% launching the 3D reconstruction.

close all; clear;

function 11_Comb_im_before_reconstr

Data = load('Cleaned_results_HT.txt');
Grain_ID = unique(Data(:,7));

% Save in separate folders the images relative to different grains
rmdir('Comb_blobs_final_grouped', 's');
mkdir('Comb_blobs_final_grouped');
for bb = 1:size(Grain_ID, 1)
   name_dir = sprintf('Comb_blobs_final_grouped/Grain_%02i', bb);
   mkdir(name_dir);
end

for ii = 1:size(Grain_ID, 1)
    [r,] = find(Data(:,7) == Grain_ID(ii));
    Data_grain = zeros(size(r,1), size(Data,2));
    Data_grain(:,:) = Data(r,:);
    % Find the unique (omega, lambda) combinations
    OL = unique(Data_grain(:, 5:6), 'rows');
    OL_flagged = zeros(size(OL,1), size(OL,2) + 1);
    OL_flagged(:,1:size(OL,2)) = OL(:,:);
    Omega = unique(OL(:,1));
    for jj = 1:size(Omega, 1)
        [r_angle, ] = find(OL(:,1) == Omega(jj));
        % If more than one image per projection, combine them
        if size(r_angle, 1) > 1
            OL_flagged(r_angle, 3) = 1;
            Image_IM = OL(r_angle,2);
            IM_comb = zeros(96, 96);
            IM_comb_final = zeros(96, 96);
            % Load images, combine them, save with name of the first in a
            % separate folder
            for kk = 1:size(r_angle, 1)
                IM_name = sprintf('Comb_blobs_final/Comb_final_%03i_%03i.png', Omega(jj), Image_IM(kk));
                IM = imread(IM_name);
                IM_scaled = IM./max(max(IM));
                for mm = 1:96
                    for nn = 1:96
                        IM_comb(mm,nn) = IM_comb(mm,nn) + IM_scaled(mm,nn);
                    end
                end
                [r_el, ] = find(OL(:,1) == Omega(jj) & OL(:,2) == Image_IM(kk));
            end
            % In the final shape, pixels have value one if they are present
            % in at least half of the components
            for mm = 1:96
                for nn = 1:96
                    if IM_comb(mm,nn) >= (size(r_angle, 1))*0.5;
                        IM_comb_final(mm,nn) = 1;
                    end
                end
            end
            % Fill holes and get rid of disconnected pixels
            [biggest_blob] = 9_Funct_clean_blobs(IM_comb_final);
            % Save the combined image. As an image name, use the ID of the
            % first image
            output = sprintf('Comb_blobs_final_grouped/Grain_%02i/Comb_final_%03i_%02i.png', ii, Omega(jj), Image_IM(1));
            if size(biggest_blob,1) > 0
                imwrite(biggest_blob, output);
            end
        end
    end
    % Copy the images that haven't been combined
    [r_uncombined, ] = find(OL_flagged(:,3) == 0);
    for aa = 1:size(r_uncombined, 1)
        Proj_n = OL_flagged(r_uncombined(aa), 1);
        img_ID = OL_flagged(r_uncombined(aa), 2);
        source_file = sprintf('Comb_blobs_final/Comb_final_%03i_%03i.png', Proj_n, img_ID);
        destination = sprintf('Comb_blobs_final_grouped/Grain_%02i/Comb_final_%03i_%03i.png', ii, Proj_n, img_ID);
        copyfile(source_file, destination);
    end
end

% Save data relative to the various images, using the format required by the
% 3D reconstruction code: Column 1: Projection #; C2: image num; C3: X_cm;
% C4: Y_cm; C5: tag

num_img = 0;
for ij = 1:size(Grain_ID, 1)
    imagefiles = dir(sprintf('Comb_blobs_final_grouped/Grain_%02i/', ij));
    nfiles = length(imagefiles);    % Number of files found
    path_dir = sprintf('Comb_blobs_final_grouped/Grain_%02i/', ij);

    for jk = 1:nfiles-2
        currentfilename = imagefiles(jk+2).name;
        path_image = strcat(path_dir, currentfilename);
        % Get image characteristics from the file name
        Im_data = sscanf(currentfilename,'%*4c_%*5c_%u_%u.png', 2);
        Angle = Im_data(1);
        IM_id = Im_data(2);
        A = imfill(double(imread(path_image)), 'holes');
        IM = int8(A);
        IM_scaled = IM./max(max(IM));
        Area = nnz(IM_scaled);
        centr = regionprops(IM_scaled,'centroid');
        centroids = cat(1, centr.Centroid);
        if isempty(centroids) == 0
            X_CM = centroids(1);
            Y_CM = centroids(2);
            num_img = num_img + 1;
            Data_img(num_img, 1:5) = [Angle, IM_id, X_CM, Y_CM, ij];
        end
    end
end

% Show, for each grain, how many projections we have
figure; histogram(Data_img(:,5));
xlabel('Grain number'); ylabel('Images available');
title('Number of images per grain');

% Write all data to a file
dlmwrite('3D_reconstr_input.txt', Data_img, 'delimiter',' ', 'precision', 6);
