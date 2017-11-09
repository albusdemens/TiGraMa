% Alberto Cereser, 20 October 2015
% Technical University of Denmark, alcer@fysik.dtu.dk

% Script to group results from the HT and get a flavour of the
% corresponding reconstruction

clear; close all;

% Dimensions of the ellipsoid used to group the data
a = 6;  % X coord
b = 6;  % Y coord
c = 6; % Z coord
% Even if the points are distributed forming ellipsoids, from the 3D
% reconstruction we have that it's better to use large parameters to avoid
% redundancies (same grain with different tags), at least when working
% with images with more than 1000 pixels

11_Plot_and_group_CM_Murofi(a, b, c);
11_Clean_result_HT(a, b, c);
11_Comb_im_before_reconstr;
%[tag] = Voxels_tagging_final;
%figure; mesh(tag);
