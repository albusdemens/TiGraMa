% Alberto Cereser, 3 January 2016
% Technical University of Denmark, alcer@fysik.dtu.dk

% This script returns the properties of the images (isolated blobs) in a
% folder. Process done for all the area ranges

close all; clear;

directory_name = sprintf('Isolated_blobs_1000_all');
files = dir(directory_name);
fileIndex = find(~[files.isdir]);
image_values = zeros(length(fileIndex), 4);
line_num = 0;
for i = 1:length(fileIndex)
    line_num = line_num + 1;
    filename = files(fileIndex(i)).name;
    val = sscanf(filename, '%*8c_%*4c_%u_%u_%u.png');
    image_values(line_num, 1:3) = val(:);
    image_values(line_num, 4) = 1;
end

disp('1 Done');
dlmwrite('Properties_isolated_blobs_1000.txt', image_values, 'delimiter', '\t');
clear;

directory_name = sprintf('Isolated_blobs_100_all');
files = dir(directory_name);
fileIndex = find(~[files.isdir]);
image_values = zeros(length(fileIndex), 4);
line_num = 0;
for i = 1:length(fileIndex)
    line_num = line_num + 1;
    filename = files(fileIndex(i)).name;
    val = sscanf(filename, '%*8c_%*4c_%u_%u_%u.png');
    image_values(line_num, 1:3) = val(:);
    image_values(line_num, 4) = 2;
end

disp('2 Done');
dlmwrite('Properties_isolated_blobs_100.txt', image_values, 'delimiter', '\t');
clear;

directory_name = sprintf('Isolated_blobs_100-1000_all');
files = dir(directory_name);
fileIndex = find(~[files.isdir]);
image_values = zeros(length(fileIndex), 4);
line_num = 0;
for i = 1:length(fileIndex)
    line_num = line_num + 1;
    filename = files(fileIndex(i)).name;
    val = sscanf(filename, '%*8c_%*4c_%u_%u_%u.png');
    image_values(line_num, 1:3) = val(:);
    image_values(line_num, 4) = 3;
end

disp('3 Done');
dlmwrite('Properties_isolated_blobs_100-1000.txt', image_values, 'delimiter', '\t');
clear;

directory_name = sprintf('Isolated_blobs_100-500_all');
files = dir(directory_name);
fileIndex = find(~[files.isdir]);
image_values = zeros(length(fileIndex), 4);
line_num = 0;
for i = 1:length(fileIndex)
    line_num = line_num + 1;
    filename = files(fileIndex(i)).name;
    val = sscanf(filename, '%*8c_%*4c_%u_%u_%u.png');
    image_values(line_num, 1:3) = val(:);
    image_values(line_num, 4) = 4;
end

disp('4 Done');
dlmwrite('Properties_isolated_blobs_100-500.txt', image_values, 'delimiter', '\t');
