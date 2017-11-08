% Alberto Cereser, 30 Sep 2014
% Department of Physics, Technical University of Denmark
% alcer@fysik.dtu.dk

% This code calculates the OB for each projection starting from those collected
% At the beginning and at the end.
clear;

number_projections = 61;
imgs_per_proj = 2422; % Number of images per projection

% Lets load the OB images relative to the various wavelengths

OB_start = zeros(imgs_per_proj, 512, 512);
OB_end = zeros(imgs_per_proj, 512, 512);

l = 0;
for k = 1:imgs_per_proj
	  s = sprintf('/data/alcer/Beamtimes/Data_BL18_Feb2015/MCP_detector/13_OB_Mar06/Corrected/13_OB_000_%05i.fits',k-1);
    t = sprintf('/data/alcer/Beamtimes/Data_BL18_Feb2015/MCP_detector/16_OB_Mar08/Corrected/16_OB_000_%05i.fits',k-1);
    l = l + 1;
    OB_start(l,:,:) = fitsread(s);
    OB_end(l,:,:) = fitsread(t);
end

% Then for each projection we simulate an OB stack

%OB_correction_1 = reshape(OB_correction, [1,1, size(OB_correction)]);
%OB_summed_start_1 = reshape(OB_summed_start, [1,1, size(OB_summed_start)]);

m = 0;
for j =1:number_projections
    for k = 1:imgs_per_proj
        R = zeros([512 512]);
        m = m + 1;
        OB_calculated(k,:,:) = ((squeeze(OB_start(k,:,:)) - squeeze(OB_end(k,:,:)))/number_projections)*(number_projections - (j - 1)) + squeeze(OB_end(k,:,:));
        t = sprintf('/data/alcer/Data_analysis/BL18_Mar2015/Sim_OB/OB_%03i_%05i.fits',j-1, k-1);
        R = squeeze(OB_calculated(k,:,:));
        fitswrite(R, t);
        if mod(k,20)==0
            clear R;
        end
        clear OB_calculated, t;
    end
end
