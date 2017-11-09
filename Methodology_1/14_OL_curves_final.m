% Alberto Cereser, 11 May 2016
% Technical University of Denmark, alcer@fysik.dtu.dk

% Script to unify, for the grains in the final reconstruction, the OL curves
% from the partial reconstructions. Returns a file listing the OL values
% before and after using the CM criterium

%close all; clear;

% Load the various volumes
[V1_refined, V, Vol_1000_plus, Vol_100_2000, Vol_100_1000, Vol_100_500] = 14_OL_curves_fnct;

V = V1_refined;

% OL before CM check
ID_V = unique(V);
L1_clean = zeros(1,3);
return;
for ii = 2:size(ID_V)

    %%%%% Volume 100-500 %%%%%

    jj = 0; % Counts the number of combinations ID partial reconstruction - final one
    G_ID = ID_V(ii);
    for aa = 1:size(V,1)
        for bb = 1:size(V,2)
            for cc = 1:size(V,3)
                if V(aa, bb, cc) == G_ID
                    if Vol_100_2000(aa, bb, cc) > 0
                        G_100_2000 = Vol_100_2000(aa, bb, cc);
                        jj = jj + 1;
                        % For each voxel, list the relative grain in V and in
                        % the partial reconstruction
                        L1(jj,1) = 1;
                        L1(jj,2) = G_ID;
                        L1(jj,3) = G_100_2000;
                    end
                end
            end
        end
    end
    if exist('L1') == 1
        L1_partial_clean = unique(L1, 'rows');
        L1_clean = [L1_clean', L1_partial_clean']';
        clear L1_partial_clean;
        clear L1;
    end

    %%%%% Volume 1000+ %%%%%

    jj = 0; % Counts the number of combinations ID partial reconstruction - final one
    G_ID = ID_V(ii);
    for aa = 1:size(V,1)
        for bb = 1:size(V,2)
            for cc = 1:size(V,3)
                if V(aa, bb, cc) == G_ID
                    if Vol_1000_plus(aa, bb, cc) > 0
                        G_1000_plus = Vol_1000_plus(aa, bb, cc);
                        jj = jj + 1;
                        % For each voxel, list the relative grain in V and in
                        % the partial reconstruction
                        L1(jj,1) = 4;
                        L1(jj,2) = G_ID;
                        L1(jj,3) = G_1000_plus;
                    end
                end
            end
        end
    end
    if exist('L1') == 1
        L1_partial_clean = unique(L1, 'rows');
        L1_clean = [L1_clean', L1_partial_clean']';
        clear L1_partial_clean;
        clear L1;
    end

    %%%%% Volume 100-2000 %%%%%

    jj = 0; % Counts the number of combinations ID partial reconstruction - final one
    G_ID = ID_V(ii);
    for aa = 1:size(V,1)
        for bb = 1:size(V,2)
            for cc = 1:size(V,3)
                if V(aa, bb, cc) == G_ID
                    if Vol_100_2000(aa, bb, cc) > 0
                        G_100_2000 = Vol_100_2000(aa, bb, cc);
                        jj = jj + 1;
                        % For each voxel, list the relative grain in V and in
                        % the partial reconstruction
                        L1(jj,1) = 3;
                        L1(jj,2) = G_ID;
                        L1(jj,3) = G_100_2000;
                    end
                end
            end
        end
    end
    if exist('L1') == 1
        L1_partial_clean = unique(L1, 'rows');
        L1_clean = [L1_clean', L1_partial_clean']';
        clear L1_partial_clean;
        clear L1;
    end


    %%%%% Volume 100-1000 %%%%%

    jj = 0; % Counts the number of combinations ID partial reconstruction - final one
    G_ID = ID_V(ii);
    for aa = 1:size(V,1)
        for bb = 1:size(V,2)
            for cc = 1:size(V,3)
                if V(aa, bb, cc) == G_ID
                    if Vol_100_1000(aa, bb, cc) > 0
                        G_100_1000 = Vol_100_1000(aa, bb, cc);
                        jj = jj + 1;
                        % For each voxel, list the relative grain in V and in
                        % the partial reconstruction
                        L1(jj,1) = 2;
                        L1(jj,2) = G_ID;
                        L1(jj,3) = G_100_1000;
                    end
                end
            end
        end
    end
    if exist('L1') == 1
        L1_partial_clean = unique(L1, 'rows');
        L1_clean = [L1_clean', L1_partial_clean']';
        clear L1_partial_clean;
        clear L1;
    end
end

L1_clean_final = unique(L1_clean, 'rows');
dlmwrite('Comb_V1_refined.txt', L1_clean_final, 'delimiter', '\t');
