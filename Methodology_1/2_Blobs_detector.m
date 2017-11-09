% Alberto Cereser, August 2015
% Technical University of Denmark

% Script to run Blobs_detector for a set of wavelenghts

% Select desired range of projection number
for Omega = 1:180
    % Select desired wavelength intervals
    for Lambda = 25:2391
        2_Blobs_detector_function(Omega, Lambda);
    end
end
