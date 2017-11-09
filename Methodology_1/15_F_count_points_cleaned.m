% Alberto Cereser, 5 December 2015
% Technical University of Denmark, alcer@fysik.dtu.dk

% Function calculating the number of points actually used for the HT (take out omegas
% with lots of points and superimposed values)

function OL_clean = 15_F_count_points_cleaned(OL, cut_proj)

count = histcounts(OL(:,2), 177);
max_num_el = max(count);
cut_value = cut_proj*max_num_el;
n_clean = 0;

for jj = 58%1:3:144
    [r_jj,] = find(OL(:,2) == jj-1);
    n_jj = size(r_jj, 1);
    num_el = sum(OL(:,2) == (jj-1));
    %if num_el < cut_value && num_el > 0
    if num_el > 0
        for ii = 1:n_jj
            n_clean = n_clean + 1;
            OL_clean(n_clean,:) = OL(r_jj(ii),:);
        end
    end
end
for jj = 154:6:178
    [r_jj,] = find(OL(:,2) == jj-1);
    n_jj = size(r_jj, 1);
    num_el = sum(OL(:,2) == (jj-1));
    %if num_el < cut_value && num_el > 0
    if num_el > 0
        for ii = 1:n_jj
            n_clean = n_clean + 1;
            OL_clean(n_clean,:) = OL(r_jj(ii),:);
        end
    end
end
