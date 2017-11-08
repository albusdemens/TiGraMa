% Alberto Cereser, 22 Mar 2015
% Department of Physics, Technical University of Denmark
% alcer@fysik.dtu.dk

% This script normalizes the collected signal using a rolling median
% The function is called by Correction_background.m

function Correction_background_function(j,m)

max_number_images = 2417; % Total: 2423
offset = 20;
IM=zeros([max_number_images-5 512 512]);

for k = 6:max_number_images
    m = m + 1;
    % And now the images cleaned by OB_correction_function
    u = sprintf('/data/alcer/Data_analysis/BL18_Mar2015/OB_cleaned_51_59/Cleaned_Fe_%03i_%05i.fits',j-1, k-1);
    IM(m,:,:) = fitsread(u);
end

%We use a rolling median (20 images before and 20 after)
rolling_median_IM = zeros([512 512]);
clean_images = zeros([512 512]);

for k = offset+1:max_number_images-30
    rolling_interval_IM = squeeze(IM(k-offset:1:k+offset,:,:));
    rolling_median_IM(:,:) = squeeze(median(rolling_interval_IM,1));
    clean_images(:,:) = (squeeze(IM(k,:,:))./(squeeze(rolling_median_IM(:,:))));

    v = sprintf('/data/alcer/Data_analysis/BL18_Mar2015/Images_divided_roll_med_51_59/Im_div_roll_median_Fe_%03i_%05i.fits',j-1, k+4);
    %R = squeeze(clean_images(:,:));
    fitswrite(clean_images, v);
end
