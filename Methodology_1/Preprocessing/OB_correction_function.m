% Alberto Cereser, 12 Mar 2015
% Department of Physics, Technical University of Denmark
% alcer@fysik.dtu.dk

% This function calculates, for a single projection, the images using the scaled
% OB. The function is called by OB_correction_no_median.m

function OB_correction_function(j,m)

number_images_per_projection = 2422; % Total: 2422
offset = 5;
OB=zeros([number_images_per_projection 512 512]);
IM=zeros([number_images_per_projection 512 512]);

for k = 1:number_images_per_projection
    % Lets load the OB images
    t = sprintf('/data/alcer/Data_analysis/BL18_Mar2015/OB_51_59/OB_%03i_%05i.fits',j-1, k-1);
    %t = sprintf('/Volumes/J-PARC_2015_3/02_OB_Feb26/Corrected/02_OB_000_%05i.fits', k-1);
    m = m + 1;
    OB(m,:,:) = fitsread(t);
    % And now the corrected images
    u = sprintf('/data/alcer/Data_analysis/BL18_Mar2015/Corrected_images_last_projections/15_Fe_%03i_%05i.fits',(j-52)/2, k-1);
    IM(m,:,:) = fitsread(u);
end
%median_IM = squeeze(median(IM,1));
%median_OB = squeeze(median(OB,1));
%We use a rolling median (5 images before and 5 after)

rolling_median_IM = zeros([512 512]);
rolling_median_OB = zeros([512 512]);
clean_images = zeros([512 512]);

for k = offset+1:number_images_per_projection-offset
    rolling_interval_IM = squeeze(IM(k-offset:1:k+offset,:,:));
    rolling_interval_OB = squeeze(OB(k-offset:1:k+offset,:,:));
    rolling_median_IM(:,:) = squeeze(median(rolling_interval_IM,1));
    rolling_median_OB(:,:) = squeeze(median(rolling_interval_OB,1));
    clean_images(:,:) = (squeeze(IM(k,:,:))./(squeeze(rolling_median_OB(:,:))));
    v = sprintf('/data/alcer/Data_analysis/BL18_Mar2015/OB_cleaned_51_59/Cleaned_Fe_%03i_%05i.fits',j-1, k-1);
    %R = squeeze(clean_images(:,:));
    fitswrite(clean_images, v);

    %clear rolling_interval_IM, rolling_interval_OB, rolling_median_IM, rolling_median_OB;
    %if mod(k,20)==0
    %    clear R, clean_images, v;
    %end
end
