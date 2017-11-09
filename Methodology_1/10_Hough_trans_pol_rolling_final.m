% Alberto Cereser, 30 Sep 2015
% Department of Physics, Technical University of Denmark
% alcer@fysik.dtu.dk

% Variations of Hough_trans_pol, written to consider separately different
% regions of the XY CM distribution. Points are divided along X using a
% rolling median. Alternative approach to the one using in Chain_blobs.m

% Plotting options in Funct_HT_sec_Murofi_final.m; datafile written by
% Funct_HT_pol_rolling_final.m

clear; close all;

step_size = 7;  % Should be an odd number

% Let's start by reading in the data

Data_blobs = load('Data_blobs_final.txt');
X = Data_blobs(:,4);
Y = Data_blobs(:,3);
x_min = round(min(X));
x_max = round(max(X));
y_min = round(min(Y));
y_max = round(max(Y));

%figure; scatter(Data_blobs(:,1), X);

% Calculate the HT for the different X sections
for vv = (x_min + ((step_size-1)/2)):(x_max - ((step_size-1)/2))
    Funct_HT_pol_rolling_final(vv, step_size, Data_blobs);
end
