% 23 October 2015, Alberto Cereser
% Technical University of Denmark, alcer@fysik.dtu.dk

% To speed up Blobs_similarity_Panda_function, we divide it in functions

function [angular_diff, angular_diff_1, angular_diff_2] = Shape_comparison_function(shape_n, shape_m, Int, A, A_dil, A_er, s, numBins)

I_1 = Int(shape_n,1);
IM_1 = int8(A{1, shape_n});
I_2 = Int(shape_m,1);
IM_2 = int8(A{1, shape_m});
if I_1 > I_2
    %IM_s = IM_2;
    IM_b = IM_1;
    % This is the index of the smaller element
    index = shape_m;
else
    %IM_s = IM_1;
    IM_b = IM_2;
    index = shape_n;
end

% Find centroid of the smaller blob
%s = regionprops(IM_s,'centroid');
centroids = cat(1, s{index}.Centroid);

IM_s_dil = int8(A_dil{1, index});
IM_s_er = int8(A_er{1, index});

% Find pixels of the bigger blob inside the dilated smaller one
Points_outside_big = (IM_b > 0) & (IM_s_dil == 0);
Points_inside_big = (IM_b > 0) & (IM_s_dil > 0);

% Find pixels of the bigger blob inside the eroded smaller one
Points_inside_small = (IM_b > 0) & (IM_s_er > 0);
Points_outside_small = (IM_b == 0) & (IM_s_er > 0);

% Calculate the angular range of the pixels outside the extended big
% blob
[row,col] = find(Points_outside_big == 1);
[theta, rho] = cart2pol(row-centroids(2), col-centroids(1));
binEdges = linspace(-pi, pi, numBins + 1);
[N,whichBin] = histc(theta, binEdges);
bin_size = (2*pi)/numBins;

% Calculate the angular range of the missing pixels inside the
% eroded big blob
[row_1,col_1] = find(Points_outside_small == 1);
[theta_1, rho_1] = cart2pol(row_1-centroids(2), col_1-centroids(1));
binEdges_1 = linspace(-pi, pi, numBins + 1);
[N_1,whichBin_1] = histc(theta_1, binEdges_1);

% Calculate the angular range of the pixels inside the eroded big
% blob
[row_2,col_2] = find(Points_inside_small == 1);
[theta_2, rho_2] = cart2pol(row_2-centroids(2), col_2-centroids(1));
binEdges_2 = linspace(-pi, pi, numBins + 1);
[N_2,whichBin_2] = histc(theta_2, binEdges_2);

angular_diff = (nnz(N)*bin_size)*57.2958;   % Measure of the angular
% difference between the two shapes
angular_diff_1 = (nnz(N_1)*bin_size)*57.2958;
angular_diff_2 = (nnz(N_2)*bin_size)*57.2958;
