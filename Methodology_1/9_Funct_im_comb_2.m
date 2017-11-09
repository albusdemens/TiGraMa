% Alberto Cereser, 1 Oct 2015
% Department of Physics, Technical University of Denmark
% alcer@fysik.dtu.dk

% Function to clean, for each projection, the set of combined images before
% the HT

function 9_Funct_im_comb_2(Omega, Data, rad, fid)

line_num = 0;
% Image size
Im_x = 96;
Im_y = 96;

[r, ~] = find(Data(:,1) == Omega);
IM_id = Data(r,2);
X_point = Data(r,3);
Y_point = Data(r,4);
Area_point = Data(r,5);
Data_clean = [IM_id, X_point, Y_point, Area_point];
% A matrix to track how the image labels change from original data --> tem
% folder --> final folder
Label_Matrix(:,1) = Data_clean(:,1);

% Group CM points
num_points = size(Data_clean,1);
Grouped_points = zeros(num_points, num_points);
for aa = 1:num_points
    X = Data_clean(:, 2:3);
    Y = Data_clean(aa, 2:3);
    [idx, ~] = (rangesearch(X,Y,rad));
    A = cell2mat(idx);
    neighbours = size(A,2);
    if neighbours > 1
       max_Area = max(max(Data_clean(A, 4), Data_clean(aa,4)));
       for bb = 1:neighbours
           if Data_clean(A(bb),4) > 0.5*max_Area
               Grouped_points(aa, 1) = aa;
               Grouped_points(aa, bb+1) = A(1,bb);
           end
       end
    end
end

% Clear matrix
Grouped_points( ~any(Grouped_points,2), : ) = [];  %rows
Grouped_points( :, ~any(Grouped_points,1) ) = [];  %columns

% Assign a number to each group
num_groups = nnz(unique(Grouped_points(:,1)));
max_num_el = size(Grouped_points,2) - 1;
Groups = zeros(num_groups, max_num_el + 1);
Groups(:,1) = [1:num_groups];
Groups(:,2:max_num_el + 1) = Grouped_points(:,2:max_num_el + 1);

% Move zeros to the rightmost positions
% Sort columns directly
[~,srtcol] = sort(Groups == 0,2);
% Sorted positions
sz  = size(Groups);
pos = bsxfun(@plus, (srtcol-1)*sz(1), (1:sz(1))'); % or use sub2ind
Groups1 = Groups(pos);

% Use a mask to group the elements
Mask = zeros(num_groups, num_points);
for cc = 1:num_groups
    num_im = nnz(Groups1(cc,:)) - 1;
    for dd = 1:num_im
        Mask(cc, Groups1(cc,dd + 1)) = 1;
    end
end

Groups_clean = zeros(num_groups, max_num_el + 1);
Tags = zeros(num_groups,1);
% Unify groups
for jj = 1:num_groups
    AA = Mask(jj,:);
    [r,c] = find(AA == 1);
    num_el_AA = size(r,2);
    for kk = 1:num_groups
        if jj ~= kk
            common_el = 0;
            BB = Mask(kk,:);
            [r1,c1] = find(BB == 1);
            num_el_BB = size(r1,2);
            Min_num_el = min(num_el_AA, num_el_BB);

            % Denote what line has less and more images
            if num_el_AA == Min_num_el
                M_ind = kk;
                c_short = c;
                c_long = c1;
            else
                M_ind = jj;
                c_short = c1;
                c_long = c;
            end
            for ll = 1:Min_num_el
                [row, ~] = find(c_long == c_short(ll));
                if nnz(row) > 0
                    common_el = common_el + 1;
                end
            end
            if common_el == Min_num_el
                Tags(jj,1) = M_ind;
                Tags(kk,1) = M_ind;
            end
        end
    end
end

Unique_tags = unique(Tags(:,1));
num_unique_tags = nnz(Unique_tags);
l_n = 0;
for vv = 1:num_unique_tags
    l_n = l_n + 1;
    Groups_clean(l_n, 1) = Unique_tags(vv,1);
    [r2,~] = find(Groups1(:,1) == Unique_tags(vv,1));
    images_max = nnz(Groups1(r2,:));
    Groups_clean(l_n,2:images_max) = Groups1(r2,2:images_max);
end

% Clean matrix
final_lines = nnz(Groups_clean(:,1));
Groups_clean_cl = zeros(final_lines, max_num_el + 1);
for uu = 1:num_groups
    if Groups_clean(uu,1) > 0
        Groups_clean_cl(uu,:) = Groups_clean(uu,:);
    end
end

% Order elements and spot possible repetitions
CC = Groups_clean_cl(:, 2:max_num_el + 1);
CC = sort(CC,2);
% Sort columns directly
[~,srtcol1] = sort(CC == 0,2);
% Sorted positions
sz1  = size(CC);
pos1 = bsxfun(@plus, (srtcol1-1)*sz1(1), (1:sz1(1))'); % or use sub2ind
CC = CC(pos1);
CC_unique = unique(CC, 'rows');
nu_el = size(CC_unique, 1);
Groups_clean_1 = zeros(nu_el, max_num_el + 1);
Groups_clean_1(:,1) = [1:nu_el]';
Groups_clean_1(:,2:max_num_el + 1) = CC_unique(:,:);

% Clean zero rows
Groups_clean_2 = zeros(size(Groups_clean_1));
num_clean_2 = 0;
for aaa = 1:size(Groups_clean_1, 1)
    if Groups_clean_1(aaa, 2) > 0
        num_clean_2 = num_clean_2 + 1;
        Groups_clean_2(num_clean_2, :) = Groups_clean_1(aaa, :);
    end
end
% Remove the zero rows
Groups_clean_2( ~any(Groups_clean_2,2), : ) = [];
Groups_clean_1 = Groups_clean_2;

% For each group, substitute indices with image numbers
Groups_clean_im = zeros(size(Groups_clean_1));
Groups_clean_im(:, 1) = Groups_clean_1(:, 1);
num_groups_clean = size(Groups_clean_1,1);
for ww = 1:num_groups_clean
    el_group = nnz(Groups_clean_1(ww,:));
    for zz = 2:el_group
        Groups_clean_im(ww, zz) = Data_clean(Groups_clean_1(ww, zz), 1);
    end
end

% Keep track of the image numbers in Label_Matrix
for ww = 1:num_groups_clean
    el_group = nnz(Groups_clean_1(ww,:));
    for zz = 2:el_group
        [r_lm, ~] = find(Label_Matrix(:, 1) == Groups_clean_im(ww, zz));
        tags_already_assigned = nnz(Label_Matrix(r_lm, :));
        Label_Matrix(r_lm, tags_already_assigned + 1) = Groups_clean_im(ww, 1);
    end
end

% Since some of the original images now have multiple tags, make a new
% matrix to keep track of things
max_num_tags = size(Label_Matrix, 2);
Assigned_group_tag = unique(Label_Matrix(:,2:max_num_tags));
Label_Matrix_comb = zeros(size(Assigned_group_tag,1)-1,2);
Label_Matrix_comb(:,1) = Assigned_group_tag(2:size(Assigned_group_tag,1),1);

% For each data combinatin step, save temporary results in a matrix
combination_count_1 = 0;  % Images combination counter

% For each group, combine the relative images
for ww = 1:num_groups_clean
    line_num = line_num + 1;
    IM_comb = zeros(Im_x, Im_y);
    IM_comb_final = zeros(Im_x, Im_y);
    el_group = nnz(Groups_clean_1(ww,:));
    images_tag = Groups_clean_1(ww,2:el_group);
    for zz = 1:el_group - 1
        IM = imread(sprintf('Comb_blobs/Comb_blobs_all/Comb_blob_%03i_%03i.png', Omega, Data_clean(images_tag(zz), 1)));
        IM_scaled = IM./max(max(IM));
        for mm = 1:Im_x
            for nn = 1:Im_y
                IM_comb(mm,nn) = IM_comb(mm,nn) + IM_scaled(mm,nn);
            end
        end
    end
    for mm = 1:Im_x
        for nn = 1:Im_y
            if IM_comb(mm,nn) > (el_group - 1)*0.5;
                IM_comb_final(mm,nn) = 1;
            end
        end
    end
    % Fill holes and avoid disconnected pixels; find area and CM
    [biggest_blob, Area_comb, X_comb, Y_comb] = 9_Funct_clean_blobs(IM_comb_final);

    output = sprintf('Comb_blobs_temp/Comb_final_%03i_%03i.png', Omega, ww);
    imwrite(biggest_blob, output);
    Final_data(line_num,1:5) = [Omega, ww, X_comb, Y_comb, Area_comb];
    for zz = 1:el_group - 1
        combination_count_1 = combination_count_1 + 1;
        Image_combination_1(combination_count_1, 1:3) = [Omega, Data_clean(images_tag(zz), 1), ww];
        % Combined data stored in matrix: Omega, image num, tag num
    end
end

%%%% Final cleaning and processing %%%%
% Check the CM values and see if we can combine some
Final_data_clean = zeros(num_groups_clean, max_num_el + 1);

for aa = 1:num_groups_clean
    X = Final_data(:, 3:4);
    Y = Final_data(aa, 3:4);
    [idx, ~] = (rangesearch(X,Y,rad));
    A = cell2mat(idx);
    neighbours = size(A,2);
    if neighbours > 1
       max_Area = max(max(Final_data(A, 5), Final_data(aa,5)));
       for bb = 1:neighbours
           if Final_data(A(bb),5) > 0.7*max_Area
               Final_data_clean(aa, 1) = aa;
               Final_data_clean(aa, bb+1) = A(1,bb);
           end
       end
    end
end

IM_progr_num = 0;       % Incremental variable used to call images in their
% final version
if nnz(Final_data_clean) > 0
    % Order elements and spot possible repetitions
    CC_fin = Final_data_clean(:, 2:max_num_el + 1);
    CC_fin = sort(CC_fin,2);
    % Sort columns directly
    [~,srtcol_fin] = sort(CC_fin == 0,2);
    % Sorted positions
    sz_fin  = size(CC_fin);
    pos_fin = bsxfun(@plus, (srtcol_fin-1)*sz_fin(1), (1:sz_fin(1))'); % or use sub2ind
    CC_fin = CC_fin(pos_fin);
    CC_unique_fin = unique(CC_fin, 'rows');
    nu_el_fin = size(CC_unique_fin, 1);
    Groups_clean_fin = zeros(nu_el_fin, max_num_el + 1);
    Groups_clean_fin(:,1) = [1:nu_el_fin]';
    Groups_clean_fin(:,2:max_num_el + 1) = CC_unique_fin(:,:);

    % Combine the similar elements. First, reduce Groups_clean_fin to the nnz
    % rows
    Groups_clean_fin( ~any(Groups_clean_fin,2), : ) = [];  %rows
    Groups_clean_fin( :, ~any(Groups_clean_fin,1) ) = [];  %columns

    % Find the rows with image IDs
    [r_fin, ~] = find(Groups_clean_fin(:,2) > 0);
    num_l_fin = size(r_fin, 1);
    Data_final_comb = zeros(num_l_fin, 5);   % Matrix to store data on the combined blobs
    for gg = 1:num_l_fin
        num_im_line = nnz(Groups_clean_fin(r_fin(gg),:));

        IM_comb_tot = zeros(Im_x, Im_y);
        IM_comb_tot_final = zeros(Im_x, Im_y);
        % Combine the images
        for hh = 2:num_im_line
            img_num = Groups_clean_fin(r_fin(gg), hh);
            IM_to_comb = imread(sprintf('Comb_blobs_temp/Comb_final_%03i_%03i.png', Omega, img_num));
            IM_to_comb = IM_to_comb./max(max(IM_to_comb));
            for mm = 1:Im_x
                for nn = 1:Im_y
                    IM_comb_tot(mm,nn) = IM_to_comb(mm,nn) + IM_comb_tot(mm,nn);
                end
            end
        end
        for mm = 1:Im_x
            for nn = 1:Im_y
                if IM_comb_tot(mm,nn) > (num_im_line - 1)*0.5;
                    IM_comb_tot_final(mm,nn) = 1;
                end
            end
        end
        % Fill holes and avoid disconnected pixels; find CM and area
        if nnz(IM_comb_tot_final) > 0
            [IM_final_clean, Area_fin, X_fin, Y_fin] = 9_Funct_clean_blobs(IM_comb_tot_final);
            Data_final_comb(gg, :) = [Omega, gg, X_fin, Y_fin, Area_fin];
            % Tag combined images in Label_matrix_comb
            for hh = 2:num_im_line
                img_num = Groups_clean_fin(r_fin(gg), hh);
                [r_lm_comb, c_lm_comb] = find(Label_Matrix_comb(:,1) == img_num);
                Label_Matrix_comb(r_lm_comb, 2) = gg;
            end

            IM_progr_num = IM_progr_num + 1;
            output_comb = sprintf('Comb_blobs_final/Comb_final_%03i_%03i.png', Omega, IM_progr_num);
            imwrite(IM_final_clean, output_comb);
        end
    end
    [~, y_fin] = size(Groups_clean_fin);
    Combined_im = nonzeros(unique(Groups_clean_fin(:, 2:y_fin)));

    num_untagged = 0;
    for ab = 1:size(Label_Matrix_comb,1)
        if Label_Matrix_comb(ab, 2) == 0
            num_untagged = num_untagged + 1;
            Untagged_im(num_untagged, 1) = Label_Matrix_comb(ab,1);
        end
    end

    % Keep track of the various combining steps
    Image_combination_2 = zeros(size(Image_combination_1,1), size(Image_combination_1,2) + 1);
    Image_combination_2(:,1:3) = Image_combination_1(:,1:3);
    Tag_Num = 0;
    for ab = 1:size(Groups_clean_fin,1)
        if Groups_clean_fin(ab,2) > 0
            Tag_Num = Tag_Num + 1;
            Num_Tagged_El = nnz(Groups_clean_fin(ab,:));
            for bc = 2:Num_Tagged_El
                [Tagged_row, ~] = find(Image_combination_1(:,3) == Groups_clean_fin(ab,bc));
                Image_combination_2(Tagged_row, 4) = Tag_Num;
            end
        end
    end

    % Copy to the final folder the images that haven't been combined in the
    % last step. In the loop, consider only these elements
    for kl = 1:size(Untagged_im,1)
        img = Untagged_im(kl,1);
        img_indx = find(Untagged_im(:,1) == img);
        IM_progr_num = IM_progr_num + 1;

        source_file = sprintf('Comb_blobs_temp/Comb_final_%03i_%03i.png', Omega, img_indx);
        destination = sprintf('Comb_blobs_final/Comb_final_%03i_%03i.png', Omega, IM_progr_num);
        copyfile(source_file, destination);
    end

    % Include information on these images in Image_combination_2. To do so
    % copy, rescaled, the tags in Image_combination_2 that haven't been
    % combined
    [r_untagged, ~] = find(Image_combination_2(:,4) == 0);
    Last_tag_used = max(Image_combination_2(:,4));
    for cd = 1:size(r_untagged, 1)
        Gr_num = unique(Image_combination_2(r_untagged, 3));
        for de = 1:size(Gr_num, 1)
            if Image_combination_2(r_untagged(cd),3) == Gr_num(de)
                Image_combination_2(r_untagged(cd),4) = de + Last_tag_used;
            end
        end
    end
else    % In case we don't have images to combine in the second step
    % Keep track of the combined images
    Image_combination_2 = zeros(size(Image_combination_1));
    Image_combination_2(:, 1:3) = Image_combination_1(:, 1:3);
    Image_combination_2(:, 4) = Image_combination_1(:, 3);
    for mn = 1:size(Label_Matrix_comb,1)
        if Label_Matrix_comb(mn, 2) == 0
            Initial_ID = mn;
            IM_progr_num = IM_progr_num + 1;
            source_file = sprintf('Comb_blobs_temp/Comb_final_%03i_%03i.png', Omega, Initial_ID);
            destination = sprintf('Comb_blobs_final/Comb_final_%03i_%03i.png', Omega, IM_progr_num);
            copyfile(source_file, destination);
        end
    end
end

% Find the images that haven't been combined in the first step; save their
% info and the relative images
Num_images_prelim = size(Image_combination_2, 1);
for lm = 1:size(Label_Matrix, 1)
    if Label_Matrix(lm, 2) == 0
        Initial_ID = Label_Matrix(lm,1);
        IM_progr_num = IM_progr_num + 1;
        % Store data in matrix
        Num_images_prelim = Num_images_prelim + 1;
        Image_combination_2(Num_images_prelim, 1) = Omega;
        Image_combination_2(Num_images_prelim, 2) = Initial_ID;
        Image_combination_2(Num_images_prelim, 3) = IM_progr_num;
        Image_combination_2(Num_images_prelim, 4) = IM_progr_num;
        % Save images
        source_file = sprintf('Comb_blobs/Comb_blobs_all/Comb_blob_%03i_%03i.png', Omega, Initial_ID);
        destination = sprintf('Comb_blobs_final/Comb_final_%03i_%03i.png', Omega, IM_progr_num);
        copyfile(source_file, destination);
    end
end

% Write the lambda values used to make each image
fprintf(fid, '%2d %2d %2d %2d\n', Image_combination_2');
