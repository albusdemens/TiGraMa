% Alberto Cereser, 10 September 2015
% Technical University of Denmark, alcer@fysik.dtu.dk

% Script to group blobs, using the list of lambda values from list_minima_1000_NN.txt. The grouping is done using angular parameters

close; clear all;

fid = fopen('Lambda_IM_before_HT.txt', 'wt');
list_n = 0;
for Omega = 1:49
    [list_n] = Track_lambda_Panda_function(Omega, list_n, fid);
end
fclose(fid);

% To divide in separate folders images relative to different proj, use
% commands_sim_diff_lambda.sh
