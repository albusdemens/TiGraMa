% Alberto Cereser, 22 Mar 2015
% Department of Physics, Technical University of Denmark
% alcer@fysik.dtu.dk

% This script normalizes the collected signal using a rolling median
clear;
number_projections = 60; 
m = 0;

for j = 1:number_projections
    Correction_background_function(j,m)
end
