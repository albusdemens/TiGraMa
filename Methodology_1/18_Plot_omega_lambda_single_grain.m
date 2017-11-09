% Alberto Cereser, 25 October 2015
% Technical Univerity of Denmark, alcer@fysik.dtu.dk

% Script to plot the omega, lambda distribution for a single grain

function [Umax2, y1_all] = 18_Plot_omega_lambda_single_grain(Grain_num)

l_shift = 0.0731;   % Number from Study_lambda_shift.m, where the distance
                    % Expected peaks - exp values is calculated
cut = 0.3;          % The only projections considered are those with # el <
                    % max*cut

OL = load('OL_final_OK.txt');
% Calculate the average number of combinations per grain
%Avg = mean(histcounts(OL(:,1), 22));
% Start by plotting the distribution of combinations per grain
%figure; histogram(OL(:,1));
%title('Number of (omega, lambda) combinations per grain (big ones)');
%xlabel('Grain ID');
%ylabel('Number of combinations'); hold on;
%plot(xlim,[Avg,Avg],'r--','LineWidth',2);

% Plot omega, lambda for the considered grain
[r, ] = find((OL(:,1) == Grain_num) & (OL(:,4) > 1.6)); % Second condition to avoid
% ultra crowded region
Omega_grain = OL(r,2);
Lambda_grain = OL(r,4);

%figure; scatter(Omega_grain, Lambda_grain);
%title('Omega, lambda distribution for a single grain');
%xlabel('Omega'); ylabel('Lambda');

% Return the distribution of points along omega
figure; histogram(Omega_grain, 144);
str = sprintf('Number of extinction spots per projection, grain %i', Grain_num);
title(str);
xlabel('Projection number'); ylabel('Number of extinction spots');
% Temporary get rid of the projections with too many points
%count = histcounts(Omega_grain, 144);
%max_num_el = max(count);
%cut_value = cut*max_num_el;
%OL_clean = zeros(size(Omega_grain, 1), 2);

%for ii = 1:size(Omega_grain, 1)
%    num_el = sum(Omega_grain(:) == Omega_grain(ii));
%    if num_el < cut_value
%        OL_clean(ii,1) = Omega_grain(ii);
%        OL_clean(ii,2) = Lambda_grain(ii);
%    end
%end
OL_clean(:,1) = Omega_grain(:);
%OL_clean(:,2) = Lambda_grain(:) + l_shift;

%OL_clean(:,2) = Lambda_grain(:)*1.025+0.0731/2.;
OL_clean(:,2) = Lambda_grain(:)*1.015+l_shift*0.8; %reasonable

% Find the grain orientation from the (omega, lambda) values
Umax = 0;
rmax = 0;
[Nmax,Umax1,rmax1,~,~]=indexToF(OL_clean(:,1),OL_clean(:,2),1,Umax,rmax);
% Plot results of the first step of the fit
[~,~,~,~,~]=indexToF(OL_clean(:,1),OL_clean(:,2),2,Umax1,rmax1);

%rmax1 = [-0.0142   -0.2642   -0.0142];
%rmax1=[0 0 0];
Umax1 =r2U(rmax1);
[Nmax,Umax2,rmax2,~,~]=indexToF(OL_clean(:,1),OL_clean(:,2),3,Umax1,rmax1);
% Plot results of the second step of the fit
[Nmax,~,rmax3,omega_calc,lambda_calc]=indexToF(OL_clean(:,1),OL_clean(:,2),2,Umax2,rmax2);

p = polyfit(omega_calc,lambda_calc,7);
% Consider only l > 2 AA; group points around the fitting curves
n_2 = 0;
for aa = 1:size(OL_clean, 1)
    %if OL_clean(aa, 2) > 2
    n_2 = n_2 + 1;
    OL_clean_2(n_2, :) = OL_clean(aa,:);
    %end
end
% figure; scatter(OL_clean_2(:,1), OL_clean_2(:,2) + l_shift, 'b*'); hold on;
% scatter(omega_calc, lambda_calc, 'r.');
% xlabel('Omega'); ylabel('Lambda'); title('Experimental data and expected values');

% Detect the different curves
OL_c = zeros(size(omega_calc,1), 3);
OL_c(:,1) = omega_calc(:,1);
OL_c(:,2) = lambda_calc(:,1);
curve_n = 1;
for ii = 2:size(omega_calc, 1)
    if (omega_calc(ii - 1,1) < omega_calc(ii,1)) && (omega_calc(ii,1) - omega_calc(ii-1,1) == 3)
        OL_c(ii,3) = curve_n;
    else
        curve_n = curve_n + 1;
        OL_c(ii,3) = curve_n;
    end
end

% Find the equation of each calculated curve
x1 = linspace(0,180,180);
y1_all = zeros(curve_n, 180);
for jj = 1:curve_n
    [r_jj, ] = find(OL_c(:,3) == jj);
    X = OL_c(r_jj,1);
    Y = OL_c(r_jj,2);
    p = polyfit(X,Y,2);
    y1 = polyval(p, x1);
    y1_all(jj, :) = y1(1,:);
end

figure; axis([0 180 2 4.2]); hold on;
scatter(OL_clean_2(:,1), OL_clean_2(:,2), 'b*'); hold on;
%scatter(OL_c(:,1), OL_c(:,2), 'b*'); hold on;
%scatter(omega_calc, lambda_calc, 'r.');
for jj = 1:curve_n
    plot(x1, y1_all(jj, :)); hold on;
    boundedline(x1, y1_all(jj, :), 0.1, 'alpha'); hold on;
end
xlabel('Omega'); ylabel('Lambda'); title('Experimental data and fitting G vectors');

% Find the experimental points that are close to the curves (inside the band)
n_tag = 0;
for i = 1:curve_n
    for j = 1:size(y1_all, 2)
        for k = 1:size(OL_clean_2, 1)
            dist = abs(y1_all(i, j) - (OL_clean_2(k,2)));
            if abs(j - (OL_clean_2(k,1))) < 2 && dist < 0.1
                disp(dist);
                n_tag = n_tag + 1;
                disp(n_tag);
                Ol_tag(n_tag, 1) = OL_clean_2(k,1);
                Ol_tag(n_tag, 2) = OL_clean_2(k,2);
                Ol_tag(n_tag, 3) = i;
                Ol_tag(n_tag, 4) = dist;
            end
        end
    end
end

figure; axis([0 180 2 4.2]); hold on;
for jj = 1:curve_n
    plot(x1, y1_all(jj, :)); hold on;
    boundedline(x1, y1_all(jj, :), 0.1, 'alpha'); hold on;
end
scatter(Ol_tag(:,1), Ol_tag(:,2), 'r*');

% Save the important solutions in a single plot
% mkdir('Omega_lambda_plots');
% h = figure('visible','off');
% subplot(2,2,1);
% scatter(OL_clean(:,1), OL_clean(:,2), '*');
% title('Omega, lambda for a single grain');
% xlabel('Omega'); ylabel('Lambda');
% hold on;
% subplot(2,2,2);
% scatter(OL_clean(:,1), OL_clean(:,2), '*');
% title('Omega, lambda for a single grain');
% xlabel('Omega'); ylabel('Lambda');
% hold on;
% scatter(omega_calc, lambda_calc, '.');
% subplot(2,2,3);
% axis([0 180 2 4.2]); hold on;
% scatter(OL_clean_2(:,1), OL_clean_2(:,2), 'b*'); hold on;
% %scatter(OL_c(:,1), OL_c(:,2), 'b*'); hold on;
% %scatter(omega_calc, lambda_calc, 'r.');
% for jj = 1:curve_n
%     plot(x1, y1_all(jj, :)); hold on;
%     boundedline(x1, y1_all(jj, :), 0.1, 'alpha'); hold on;
% end
% xlabel('Omega'); ylabel('Lambda'); title('Experimental and expected values');
% subplot(2,2,4);
% axis([0 180 2 4.2]); hold on;
% for jj = 1:curve_n
%     plot(x1, y1_all(jj, :)); hold on;
%     boundedline(x1, y1_all(jj, :), 0.1, 'alpha'); hold on;
% end
% scatter(Ol_tag(:,1), Ol_tag(:,2), 'r*');
% title('Points within bands');
% name = sprintf('Omega_lambda_plots/Exp_calc_grain%02i', Grain_num);
% saveas(h,name,'png');
% close all;
