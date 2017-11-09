% Alberto Cereser, 17 Nov 2015
% Technical University of Denmark, alcer@fysik.dtu.dk

% Function to calculate, for a distribution of points, the accumulator
% matrix of the HT. Called by HT_orientation.m

function [Counter_matrix_one_blob_sum] = 15_Funct_acc_matrix(OL, n_large_pixels, steps_alpha, steps_size)

x_min = round(min(OL(:,2)));    % Omega values
x_max = round(max(OL(:,2)));
y_min = round(min(OL(:,4)));    % Lambda values
y_max = round(max(OL(:,4)));

Counter_matrix = zeros(n_large_pixels, n_large_pixels); % We use this matrix to locate the maxima
% of the Hough transform
Counter_matrix_bin = zeros(n_large_pixels, n_large_pixels); % This is the binarized version of
% Counter_matrix, made using the threshold "cut"
Counter_matrix_one_blob_sum = zeros(n_large_pixels, n_large_pixels);

for pp = 1:size(OL(:,1))
    disp(pp/size(OL,1));
    Counter_matrix_one_blob = zeros(n_large_pixels, n_large_pixels);
    R_n = zeros(1,(steps_alpha/steps_size));
    A_n = zeros(1,(steps_alpha/steps_size));
    n = 0;
    for Alpha = (0:steps_size:steps_alpha-1)  % This is the angle of the blob CM in the sample
        den = cosd(OL(pp,2)+Alpha);
        if (abs(den)>0)  % Conditions included to avoid problems with the den
            R = abs(OL(pp,4)/den);
            n = n + 1;
            R_n(n) = R;
            A_n(n) = Alpha;
        end
    end
    len_n_one_alpha = size(A_n,2);
    M_n_one_alpha = zeros(n_large_pixels,n_large_pixels);
    for ij = 1:n
        ix = floor(A_n(ij)/360*n_large_pixels)+1;
        iy = floor(R_n(ij)/4.0538*n_large_pixels)+1;
        M_n_one_alpha(ix,iy) = 1;
    end

    for ii = 1:n_large_pixels
        for jj = 1:n_large_pixels
            Counter_matrix_one_blob_sum(ii,jj) = Counter_matrix_one_blob_sum(ii,jj) + M_n_one_alpha(ii,jj);
        end
    end
    clear M_n_one_alpha;
end
