% Alberto Cereser, 30 November 2015
% Technical University of Denmark, alcer@fysik.dtu.dk

% Script to save the reconstruction output as a vtk file, to be
% opened in ParaView

close all; clear;

Files = dir(sprintf('Grain_reconstruction/'));
n_files = length(Files);    % Number of layers
path_files = sprintf('Grain_reconstruction/');

for aa = 1:n_files-3
    currentfilename = Files(aa+3).name;
    path_image = strcat(path_files, currentfilename);
    A = double(fitsread(path_image));
    Vol(:,:,aa) = A(:,:);
end
clear A;

n_clean = size(unique(Vol));
list_grains = unique(Vol);

output = sprintf('Reconstruction.vtk');
savevtk(Vol, output);
