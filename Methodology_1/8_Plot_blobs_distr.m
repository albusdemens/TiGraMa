% Alberto Cereser, 25 October 2015
% Technical Univeristy of Denmark, alcer@fysik.dtu.dk

% Script to plot the distribution, among projections, of the extinction
% spots combined using 5_Track_lambda.m

clear; close all;

Spots = load('Data_blobs.txt');

figure;
histogram(Spots(:,1),60);
axis([0 181 -inf inf]);
title('# of extinction spots per proj');
xlabel('Projection number');
ylabel('Number of spots');

figure; scatter3(Spots(:,3), Spots(:,4), Spots(:,1)); hold on;
xlabel('X_CM'), ylabel('Y_CM'), zlabel('Z_CM');
