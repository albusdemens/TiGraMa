% Alberto Cereser, 11 may 2016
% Technical University of Denmark, alcer@fysik.dtu.dk

% Get the extiction spots list for the different grains

close all; clear;
% Load list of grains containing partial grains, generated using
% 14_OL_curves_final.m
OL = load('Comb_V1_refined.txt');

Final_ID = unique(OL(:,2));

ID_1 = load('Reconstruction_100-500/Omega_lambda_final.txt');
ID_2 = load('Reconstruction_100-1000/Omega_lambda_final.txt');
ID_3 = load('Reconstruction_100/Omega_lambda_final.txt');
ID_4 = load('Reconstruction_1000/Omega_lambda_final.txt');

ID_list = zeros(1,5);
for ii = 2:size(Final_ID,1)
    ID_fin = Final_ID(ii);
    [r,] = find(OL(:,2) == ID_fin);
    Data_partial(:,:) = OL(r,:);
    for jj = 1:size(Data_partial,1)
        Num_dataset = Data_partial(jj,1);
        ID_partial = Data_partial(jj,3);
        if Num_dataset == 1
            [r_partial,] = find(ID_1(:,1) == ID_partial);
            List_partial = zeros(size(r_partial,1), 5);
            List_partial(:,1) = ID_fin;
            List_partial(:,2:5) = ID_1(r_partial, 1:4);
            clear r_partial;
        end
        if Num_dataset == 2
            [r_partial,] = find(ID_2(:,1) == ID_partial);
            List_partial = zeros(size(r_partial,1), 5);
            List_partial(:,1) = ID_fin;
            List_partial(:,2:5) = ID_2(r_partial, 1:4);
            clear r_partial;
        end
        if Num_dataset == 3
            [r_partial,] = find(ID_3(:,1) == ID_partial);
            List_partial = zeros(size(r_partial,1), 5);
            List_partial(:,1) = ID_fin;
            List_partial(:,2:5) = ID_3(r_partial, 1:4);
            clear r_partial;
        end
        if Num_dataset == 4
            [r_partial,] = find(ID_4(:,1) == ID_partial);
            List_partial = zeros(size(r_partial,1), 5);
            List_partial(:,1) = ID_fin;
            List_partial(:,2:5) = ID_4(r_partial, 1:4);
            clear r_partial;
        end
        ID_list = [ID_list', List_partial']';
        clear List_partial;
    end
    clear Data_partial;
end

ID_list_final = unique(ID_list, 'rows');
dlmwrite('OL_comb_after_CM.txt', ID_list_final, 'delimiter', '\t');
