% Alberto Cereser, 25 October 2015
% Technical Univerity of Denmark, alcer@fysik.dtu.dk

% Script to save the omega, lambda distribution and fit for all the grains
% in a reconstruction

clear; close all;

grain_IDs = 3;%[3 21 36];    % Grains with nice curves in the OL distribution
% [3 4 8 20 21 25 29 30 33 36 48 63 153]; % These are all the grains recorded in
                                          % The 3DND slice
for i = 1:size(grain_IDs,2)
    grain = grain_IDs(i);
    [U_final, y1_all] = 18_Plot_omega_lambda_single_grain(grain);
    %output = sprintf('/Users/Alberto/Documents/Data_analysis/J-PARC_2015/Data_analysis_final/Indexing/U_matrices/U_grain_%i.txt', grain);
    %dlmwrite(output, U_final, 'delimiter',' ');
end
