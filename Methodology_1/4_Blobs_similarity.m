% Alberto Cereser, 18 September 2015
% Technical University of Denmark, alcer@fysik.dtu.dk
% Script used to run Blobs_similarity_Panda_function.m
% For different projections

clear; close all;

% Number of pixels used to enlarge the small blob and calculate the angular
% difference with the smaller one calculated suing the formula
% frame_size = sqrt(Area_small/3)*percent_size (L. 71)
percent_size = 0.2;

for Omega = 1:49
    4_Blobs_similarity_function(Omega, percent_size);
end
