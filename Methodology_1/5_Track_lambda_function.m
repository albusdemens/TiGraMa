% Alberto Cereser, 10 September 2015
% Technical University of Denmark, alcer@fysik.dtu.dk

% Function to group blobs, using the list of lambda values from list_minima_1000_NN.txt
% where we have parameters of blobs compoared using angular parameters

function [list_n] = 5_Track_lambda_function(Omega, list_n, fid)

Omega = Omega - 1;

filename = sprintf('Results_blobs_comp/Results/List_minima/list_minima_%02i.txt', Omega);
delimiterIn = ' ';
A = importdata(filename,delimiterIn);
L = size(A,1);

Lambda_1 = A(:,1);
Lambda_2 = A(:,2);
n_tot_lam_1 = size(Lambda_1,1);
Lambda_1_u = unique(Lambda_1);
n_u_lam_1 = nnz(unique(Lambda_1));
num_blobs =  max(Lambda_1);
bin = 1:num_blobs;
B = histc(Lambda_1, bin);
max_count = max(B);

% Save in a matrix the image number (first column) and the following ones
% (following columns)
List_matrix = zeros(n_u_lam_1, max_count+1);
count_unique_1 = 0;
count_tag = 0;

for ii = 1:n_u_lam_1
    indices = find(Lambda_1 == Lambda_1_u(ii));
    num_indices = size(indices, 1);
    List_matrix(ii,1) = Lambda_1_u(ii);
    List_matrix(ii,2:num_indices+1) = Lambda_2(indices);
end
unique_lambda = unique(nonzeros(List_matrix));
n_unique_lambda = size(unique_lambda, 1);

% Use a mask to group the elements
[row, col] = size(List_matrix);
Mask = zeros(n_unique_lambda, row);
for jj = 1:n_unique_lambda
    [r1, r2] = find(List_matrix == unique_lambda(jj,1));
    Mask(jj,r1) = 1;
end

% We tag the rows and then the corresponding images
List_tag_rows = zeros(row,1);
List_tag_images = zeros(n_unique_lambda, 1);

% Compare the various rows in the mask matrix. As a condition for two
% columns to be similar, we ask that they have at least half of the
% elements on the same row
tag = 0;
for kk = 1:row
    num_el = nnz(Mask(:, kk));
    [i,j]=ind2sub(size(Mask), find(Mask(:, kk)==1));
    for ll = 1:row
        if ll ~= kk
            num_el_1 = nnz(Mask(:, ll));
            [i_1,j_1]=ind2sub(size(Mask), find(Mask(:, ll)==1));
            shorter_col = min(num_el, num_el_1);
            if num_el > num_el_1
                min_num = ll;
                max_num = kk;
            else
                min_num = kk;
                max_num = ll;
            end
            inds = find(ismember(i, i_1));
            num_inds = nnz(inds);
            % Case when the two rows have the same elements
            if num_inds == shorter_col
                if List_tag_rows(kk,1) == 0 && List_tag_rows(ll,1) == 0
                    List_tag_rows(kk,1) = kk;
                    List_tag_rows(ll,1) = kk;
                elseif List_tag_rows(kk,1) > 0 && List_tag_rows(ll,1) == 0
                    assigned_tag = nnz(List_tag_rows(kk,1));
                    List_tag_rows(max_num,assigned_tag + 1) = kk;
                    List_tag_rows(min_num,1) = kk;
                elseif List_tag_rows(kk,1) > 0 && List_tag_rows(ll,1) > 0
                    assigned_tag = nnz(List_tag_rows(kk,1));
                    assigned_tag_1 = nnz(List_tag_rows(ll,1));
                    List_tag_rows(kk,assigned_tag + 1) = kk;
                    List_tag_rows(ll,1) = kk;
                end
            % Case when the two rows have different elements
            elseif num_inds > shorter_col*0.5
                if List_tag_rows(kk,1) == 0 && List_tag_rows(min_num,1) == 0
                    List_tag_rows(kk,1) = kk;
                    List_tag_rows(ll,1) = kk;
                elseif List_tag_rows(kk,1) > 0 && List_tag_rows(ll,1) == 0
                    assigned_tag = nnz(List_tag_rows(kk,1));
                    List_tag_rows(kk,assigned_tag + 1) = kk;
                    List_tag_rows(ll,1) = kk;
                elseif List_tag_rows(kk,1) > 0 && List_tag_rows(ll,1) > 0
                    assigned_tag = nnz(List_tag_rows(kk,1));
                    assigned_tag_1 = nnz(List_tag_rows(ll,1));
                    List_tag_rows(kk,assigned_tag + 1) = kk;
                    List_tag_rows(ll,assigned_tag_1 + 1) = kk;
                end
            elseif num_inds == 0
                List_tag_rows(kk,1) = kk;
            end
        end
    end
end

A = reshape(List_tag_rows,[],1);
%figure; hist(A);

% Let's now clean the mask: cancel double entries, make things more
% coincised

[ix, iy] = size(List_tag_rows);
List_tr_clean = zeros(ix,iy);
% Get rid of ripetute tags
for mm = 1:row
    AA = List_tag_rows(mm,:);
    Unique_el = nonzeros(unique(AA));
    num_unique = nnz(Unique_el);
    for nn = 1:num_unique
        List_tr_clean(mm,nn) = Unique_el(nn);
    end
end

% If a row has multiple labels, keep only the ones appearing at least twice
% in the List_tr_clean matrix
List_trc_edit = zeros(ix,iy);
for mm = 1:row
    num_el_r = nnz(List_tr_clean(mm,:));
    AAA = List_tr_clean(mm,:);
    if num_el_r > 1
        pos = 0;
        for nn = 1:num_el_r
            Occ = nnz(List_tr_clean==AAA(nn));
            if Occ > 1
                pos = pos + 1;
                List_trc_edit(mm,pos) = List_tr_clean(mm,nn);
            end
        end
    else
        List_trc_edit(mm,:) = List_tr_clean(mm,:);
    end
end

% Now we want to tag the various images, which will then be combined
Tags_im = zeros(n_unique_lambda, ix);
Tags_im(:,1) = unique_lambda;
% For each row, we tag the corresponding images
for pp = 1:ix
    num_tags_row = nnz(List_trc_edit(pp,:));
    num_images_row = nnz(List_matrix(pp,:));
    for qq = 1:num_tags_row
        for rr = 1:num_images_row
            im_num = List_matrix(pp,rr);
            [r_line, c] = find(Tags_im(:,1) == im_num);
            tags_im_num = nnz(Tags_im(r_line,:)) - 1;
            assigned_tag = List_trc_edit(pp,qq);
            if nnz(find(Tags_im(r_line,2:ix) == assigned_tag)) == 0
                Tags_im(r_line,tags_im_num + 2) = assigned_tag;
            end
        end
    end
end

% Simplify the matrix by taking out the empty column
Tags_im( :, ~any(Tags_im,1)) = [];
[ix_1, iy_1] = size(Tags_im);

% Combine and save the images with the same tag. To do so, start by finding
% the unique tags
List_tags = unique(Tags_im(:, 2:iy_1));
Lt_pos = List_tags(List_tags > 0);
num_lt_pos = nnz(Lt_pos);
maxim = max(max(Tags_im(:, 2:iy_1)));
max_count_tag = nnz(find((Tags_im(:, 2:iy_1)) == maxim));
% Matrix where we can store, for each tag, the corresponding images
Im_per_tag = zeros(num_lt_pos, max_count_tag);
Im_per_tag(:,1) = Lt_pos(:,1);

for aa = 1:num_lt_pos
    im_list = 0;
    for bb = 1:ix_1
       [num, plc] = find(Tags_im(bb, 2:iy_1) == Lt_pos(aa));
       for cc = 1:num
           im_list = im_list + 1;
           Im_per_tag(aa, im_list + 1) = Tags_im(bb,1);
       end
    end
end

% Load the file mapping image number with lambda values
file_map = sprintf('Results_blobs_comp/Results/Map_blobs/map_blobs_%02i.txt', Omega);
delimiterIn = ' ';
Map = importdata(file_map,delimiterIn);

% Combine the images listed in Im_per_tag
for dd = 1:num_lt_pos      % Separately consider each tag
    IM_comb = zeros(257,512);
    n_blobs = nnz(size(Im_per_tag(dd,:) - 1));
    for ee = 1:n_blobs
        list_n = list_n + 1;
        image_number = Im_per_tag(dd,ee + 1);
        [r_map, c_map] = find(Map(:,1) == image_number);
        lam = Map(r_map, 2);
        blob_idx = Map(r_map, 3);
        path_image = sprintf('Isolated_blobs/Isolated_blobs_%02i/Isolated_blob_%03i_%05i_%02i.png', Omega, Omega, lam, blob_idx);
        IM = double(imread(path_image));
        maxPixelValue = mean(IM(IM>0)); % find the maximum; min not working
        minPixelValue = min(min(IM)); % find the minimum
        IM = ((IM + minPixelValue)/maxPixelValue);
        IM_comb = imadd(IM_comb, IM);
        % Store the lambda values used to make each image
        List_lambda_IM(list_n, 1) = Im_per_tag(dd,1);
        List_lambda_IM(list_n, 2) = Omega;
        List_lambda_IM(list_n, 3) = lam;
    end
    [x_comb, y_comb] = size(IM_comb);
    IM_final = zeros(x_comb, y_comb);
    IM_final = imfill(double(IM_comb>0.5*n_blobs));

    [IM_L, num] = bwlabel(IM_final);
    count_pixels_per_obj = sum(bsxfun(@eq,IM_L(:),1:num));
    [~,ind] = max(count_pixels_per_obj);
    biggest_blob = double(IM_L==ind);

    % Save the final result
    output = sprintf('Comb_blobs/Comb_blob_%03i_%02i.png', Omega, Im_per_tag(dd,1));
    imwrite(biggest_blob, output);
end

% Remove zero rows in List_lambda_IM
List_lambda_IM( ~any(List_lambda_IM,2), : ) = [];

% Write the lambda values used to make each image
fprintf(fid, '%2d %2d %3d\n', List_lambda_IM');
