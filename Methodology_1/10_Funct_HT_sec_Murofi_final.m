% Alberto Cereser, 6 Jun 2015
% Department of Physics, Technical University of Denmark
% alcer@fysik.dtu.dk

% Variations of Hough_trans_pol, written to consider separately different
% regions of the XY CM distribution. Points divided along X using a rolling median

function [Data_with_CM] = 10_Funct_HT_sec_Murofi_final(Tagged_matrix_clean,...
    threshold,threshold_distance_fit,min_points, y_min, y_max, center_interval)

% Let's start by reading in the data
number_lines_txt = numel(Tagged_matrix_clean(:,1));
center = 48; % Middle of the frame

steps_alpha = 360;
steps_size = 0.01;
n_large_pixels = 100;    % This is the size of the mesh used for the Hough
                         % transform

L = 10; % This is to check that the CM of the grain isn't outside the sample
num_tags = round(max(Tagged_matrix_clean(:,4)));

% Let's count how many times we have each tag
Blobs_counter = zeros(num_tags,1);

for mm = (1:num_tags)
    counter = 0;
    for nn = (1:number_lines_txt)
        if (Tagged_matrix_clean(nn,4) == mm)
          counter = counter + 1;
        end
    end
    Blobs_counter(mm) = counter;
end

% For each section we calculate the Hough transform, which gives us the
% CM of the blob relative to the extinction spot in the sample
Data_with_CM = zeros(1,7);
for tag = 1:num_tags
    List_images = zeros(Blobs_counter(tag),6);
    count = 0;
    for ll = 1:number_lines_txt
        if Tagged_matrix_clean(ll,4) == tag
            % We make a temporary matrix where we store info about a certain
            % tag
            count = count + 1;
            List_images(count,1:6) = Tagged_matrix_clean(ll,1:6);
            % X, Y, projection, tag, number of pixels, lambda
        end
    end

    % Now we load everything we need to calculate the Hough transform
    Angle = List_images(:,3);
    %Angle = Omega*3;
    X_measured = List_images(:,1);
    Y_measured = List_images(:,2);
    Area = List_images(:,5);
    Lambda = List_images(:,6);
    Number_points = size(Angle, 1); % This is the number of points we deal with

    % Calculate the Hough transform for all points together
    ID_big = find(Area>0);
    Y_big_spots = Y_measured(ID_big);
    Angle_big_spots = Angle(ID_big);
    num_big = size(ID_big,1);
    [number_centers_big, tagged_points_big, alpha_big, R_big] = ...
        10_Funct_Hough_split_final(Y_big_spots, Angle_big_spots, ...
        center, n_large_pixels, num_big, steps_alpha, steps_size, ...
        threshold, threshold_distance_fit);

    % Store the values that are not relative to any curve
    if (size(tagged_points_big,2) == 2)
        Y_big_unfit = zeros(num_big-nnz(tagged_points_big(:,2)),1);
        Angle_big_unfit = zeros(num_big-nnz(tagged_points_big(:,2)),1);
        count_big_unfit = 0;
        for zz = 1:num_big
            if tagged_points_big(zz,2) == 0
                count_big_unfit = count_big_unfit + 1;
                Y_big_unfit(count_big_unfit) = Y_big_spots(zz);
                Angle_big_unfit(count_big_unfit) = Angle_big_spots(zz);
            end
        end
    else
        Y_big_unfit = Y_big_spots;
        Angle_big_unfit = Angle_big_spots;
    end
    % Unify all CM values
    num_CM_big = nnz(alpha_big);
    alpha_R = zeros(num_CM_big, 2);

    if num_CM_big > 0
        [min_f, max_f] = 10_Min_max_fitting_function(num_CM_big, alpha_big, R_big, center);
        for ae = 1:num_CM_big
            % Keep only fitting curves that don't go too far from the
            % considered values
            %min_big = min(nonzeros(tagged_points_big(:,2:num_CM_big+1)));
            %max_big = max(nonzeros(tagged_points_big(:,2:num_CM_big+1)));
            if ((min_f(ae) > (y_min - 5)) && (max_f(ae) < (y_max + 5)))
                % Store only the (R, alpha) values fitting at least 10% of the
                % points
                if nnz(tagged_points_big(:,ae+1)) > min_points
                    alpha_R(ae,1) = alpha_big(ae);
                    alpha_R(ae,2) = R_big(ae);
                end
            end
        end
    end

	  alpha_fit_pre = alpha_R(:,1);
	  R_fit_pre = alpha_R(:,2);
    num_nnz = 0;
    if nnz(alpha_fit_pre(:,1)) > 0
        alpha_fit = zeros(num_nnz, 1);
        R_fit = zeros(num_nnz, 1);
        for ab = 1:size(alpha_fit_pre,1)
            if alpha_fit_pre(ab,1) > 0
                num_nnz = num_nnz + 1;
                alpha_fit(num_nnz) = alpha_fit_pre(ab);
                R_fit(num_nnz) = R_fit_pre(ab);
            end
        end
    end
    num_CM_total = num_nnz;

    % For each point, we want to find what is the closest fitting line.
    % This is done in two steps: at first we check if the distance is below
    % a threshold, and then we look for the lowest distance value.
    fit_value = zeros(Number_points, num_CM_total);
    dist_curve = zeros(Number_points, num_CM_total);
    dist_curve_thresh = zeros(Number_points, num_CM_total);
    % Here we store the data after using a threshold
    dist_curve_thresh_fin = zeros(Number_points, num_CM_total);
    % Here we store the data after calculating the minimal distance
    number_candidates = zeros(Number_points);
    for aa = 1:Number_points
        for bb = 1:num_CM_total
            angular_coord = Angle(aa);
            fit_value(aa,bb) = double((R_fit(bb)*sind(alpha_fit(bb)+angular_coord)/0.055) + center);
            dist_curve(aa,bb) = abs(fit_value(aa,bb) - Y_measured(aa));
            if (dist_curve(aa,bb) < threshold_distance_fit)
                dist_curve_thresh(aa,bb) = dist_curve(aa,bb);
            end
        end
    end
    for cc = 1:Number_points
        for dd = 1:num_CM_total
            number_candidates(cc) = sum(dist_curve_thresh(cc,:)~=0,2);
            if number_candidates(cc) > 1
                index_min_dist = find(dist_curve_thresh(cc, :) > 0, 1, 'first');
                dist_curve_thresh_fin(cc,index_min_dist) = dist_curve_thresh(cc, index_min_dist);
            else
                dist_curve_thresh_fin(cc,dd) = dist_curve_thresh(cc,dd);
            end
        end
    end

    % We label each point in the Angle, Y _measured space so that we know
    % what's the closest curve
    tagged_points = zeros(Number_points,7);
    tagged_points_edit = zeros(Number_points, num_CM_total + 1);
    for ee = 1:Number_points
        for ff = 1:num_CM_total
            if dist_curve_thresh_fin(ee,ff) > 0
                tagged_points(ee,1) = Angle(ee);
                tagged_points(ee,2) = X_measured(ee);
                tagged_points(ee,3) = Y_measured(ee);
                tagged_points(ee,4) = Lambda(ee);
                tagged_points(ee,5) = ff;
                tagged_points(ee,6) = alpha_fit(ff);
                tagged_points(ee,7) = R_fit(ff);
            end
        end
        for gg = 1:num_CM_total
            if tagged_points(ee,4) == gg
                tagged_points_edit(ee,1) = tagged_points(ee,1); % Angle
                tagged_points_edit(ee,gg+1) = tagged_points(ee,3);  % Y_measured
            end
        end
    end

    % To optimize the threshold value, let's check how many empty rows we
    % have in dist_curve_thresh_fin
    percentage_labelled = nnz(dist_curve_thresh_fin)/Number_points;

    % We can now plot the points, the curve and the distance values
    figure;
    plot(Angle, Y_measured, '*'); hold on;
    f = zeros(num_CM_total, 180);
    for xx = 1:num_CM_total
        for deg = 1:180
    %        % f is the fitting function
            f(xx, deg) = double((R_fit(xx)*sind(alpha_fit(xx)+deg)/0.055) + center);
        end
        plot(1:180, f(xx,:), '.'); hold on;
        boundedline(1:180, f(xx,:), threshold_distance_fit, 'alpha');
    end
    %xlabel('Omega'); ylabel('Y detector');
    %title('Results from the HT');

    line_num = 0;
    Data_with_CM = zeros(nnz(tagged_points(:,2)),7);
    % We can now write the grains CM positions to a text file
    for lll = 1:number_lines_txt
        if size(tagged_points,1) > 0
            if tagged_points(lll,2) > 0
                line_num = line_num +1;
                Data_with_CM(line_num,1) = tagged_points(lll,1);     % Projection number
                Data_with_CM(line_num,2) = center_interval;          % Center of the considered interval
                Data_with_CM(line_num,3:5) = tagged_points(lll,2:4); % X and Y of the spot CM; lambda tagging the corresponding image
                Data_with_CM(line_num,6:7) = tagged_points(lll,6:7); % alpha and R values from HT
            end
        end
    end
end
