% Alberto Cereser, 11 Mar 2015
% Department of Physics, Technical University of Denmark
% alcer@fysik.dtu.dk

% This script corrects, for each projection, the images using the scaled OB
clear;
number_projections = 49;%61;
m = 0;

for j = 1:number_projections
    OB_correction_function(j,m)
end
