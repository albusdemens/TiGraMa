% Alberto Cereser, 18 September 2015
% Technical University of Denmark, alcer@fysik.dtu.dk
% Function to calculate, for the desired projection,
% Which blobs are similar and to group them

% Version for blobs with I in [100, 2000]

function Blobs_similarity_function(Omega, percent_size)

num_blobs = 0;

% Load all images for the considered projection
imagefiles = dir(sprintf('Isolated_blobs_%02i', Omega-1));
nfiles = length(imagefiles);    % Number of files found
numBins = 100;

for ii = 3:nfiles
   currentfilename = imagefiles(ii).name;
   path_dir = sprintf('Isolated_blobs_%02i/', Omega-1);
   path_image = strcat(path_dir, currentfilename);
   B = imfill(double(imread(path_image)), 'holes');
   % We rescale the intensity values, so that the maximum is 1
   B_scaled = B./max(max(B));
   A{ii} = B_scaled;
   % We also store the dilated and eroded version of each blob, to avoid during
   % inside the loop below. The number of pixels in the dilation/erosion process
   % depends on the size of the smaller blob
   Area_blob = nnz(A{ii});
   frame_size = round(sqrt(Area_blob/3)*percent_size);
   SE = strel('disk', frame_size);
   A_dil{ii} = imdilate(B_scaled, SE);
   A_er{ii} = imerode(B_scaled, SE);
   s{ii} = regionprops(B_scaled,'centroid');

   num_blobs = num_blobs + 1;
   out = regexp(path_image,'\d+','match');
   l_out = size(out, 2);
   lam = str2double(out(:, l_out - 1));
   idx = str2double(out(:, l_out));
   Int(num_blobs,1) = ii;           % Image number
   Int(num_blobs,2) = lam;          % Lambda
   Int(num_blobs,3) = idx;          % Blob index (there can be various blobs
                                    % for one lambda)
   Int(num_blobs,4) = Area_blob;   % Area of the blob in the considered image
end

% Write to a text file the info about all the blobs
file_blobs = fopen(sprintf('Results_blobs_comp/Results/Map_blobs/map_blobs_%02i.txt', Omega-1),'w');
fprintf(file_blobs,'%4d %4d %2d %5d\n',Int');

% Check if the two blobs have similar shape
comp_index = 0;
nfiles_true = nfiles-3;     % The first three files are not images
Comparison_list = zeros(nfiles_true*nfiles_true, 4);
for shape_n = 3:nfiles_true
        disp(shape_n);
    for shape_m = shape_n + 1:nfiles_true-2
        centroid_n = cat(1, s{shape_n}.Centroid);
        centroid_m = cat(1, s{shape_m}.Centroid);
        dist_centroids = sqrt((centroid_n(1) - centroid_m(1))^2 + (centroid_n(2) - centroid_m(2))^2);
        if dist_centroids < 50  % 50 for small blobs, 100 for the big ones

            [angular_diff, angular_diff_1, angular_diff_2] = Shape_comparison_function(shape_n, shape_m, Int, A, A_dil, A_er, s, numBins);

            if angular_diff < 45
                if angular_diff_1 < 45 && angular_diff_2 > 315
                    comp_index = comp_index + 1;
                    Comparison_list(comp_index, 1) = shape_n;
                    Comparison_list(comp_index, 2) = shape_m;
                    Comparison_list(comp_index, 3) = angular_diff;
                    Comparison_list(comp_index, 4) = angular_diff_1;
                end
            end
        end
    end
end

% Clear the matrix
Comparison_list( ~any(Comparison_list,2), : ) = [];

dlmwrite(sprintf('Results_blobs_comp/Results/List_minima/list_minima_%02i.txt', Omega-1),Comparison_list, 'delimiter',' ', 'precision', 6);
