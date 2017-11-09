% Alberto Cereser, 8 Jun 2015
% Technical University of Denmark, alcer@fysik.dtu.dk

% Function to find the min and max of the fitting sinuosoidal functions found in
% Function_HT_section_Murofi_split

function [min_f, max_f] = Min_max_fitting_function(num_CM_considered, alpha_fit, R_fit, center)

f = zeros(num_CM_considered, 180);
for xx = 1:num_CM_considered
    for deg = 1:180
        % f is the fitting function
        f(xx, deg) = double((R_fit(xx)*sind(alpha_fit(xx)+deg)/0.055) + center);
    end
end

min_f = zeros(num_CM_considered, 1);
max_f = zeros(num_CM_considered, 1);

for yy = 1:num_CM_considered
    min_f(yy) = min(f(yy,:));
    max_f(yy) = max(f(yy,:));
end