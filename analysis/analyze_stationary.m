% analyze_stationary.m
clear all;
load('EMT_ensemble_data.mat');

%% 1. Flatten the Ensemble Data
% We discard the first 20% of each run as "burn-in" equilibration time
burn_in_fraction = 0.2; 

A_flat = [];
R_flat = [];

num_runs = length(A_ensemble);
for i = 1:num_runs
    T_len = length(T_ensemble{i});
    start_idx = round(burn_in_fraction * T_len) + 1;
    
    % Extract the slice
    A_temp = A_ensemble{i}(start_idx:end);
    R_temp = R_ensemble{i}(start_idx:end);
    
    % Force them to be column vectors using (:) before stacking
    A_flat = [A_flat; A_temp(:)];
    R_flat = [R_flat; R_temp(:)];
end

%% 2. Plot the 2D Density Map (Stationary Distribution)
figure('Position', [100, 100, 700, 600]);

% Create a 2D histogram (density map)
h = histogram2(A_flat, R_flat, 'DisplayStyle', 'tile', 'ShowEmptyBins', 'on');
h.XBinLimits = [0 N];
h.YBinLimits = [0 N];
h.NumBins = [N+1 N+1]; % Bins from 0 to N

% Set the color scale to logarithmic to reveal shallow transitional basins
set(gca, 'ColorScale', 'log'); 
h.ShowEmptyBins = 'off'; % Keeps true zero values black

% Styling to make it look like a landscape heatmap
colormap(jet);
cb = colorbar;
cb.Label.String = 'Probability Density in Log scale ';
cb.Label.FontSize = 12;

xlabel('H3K4me3 (Activating / Mesenchymal)', 'FontSize', 14, 'FontWeight', 'bold');
ylabel('H3K27me3 (Repressive / Epithelial)', 'FontSize', 14, 'FontWeight', 'bold');
title('EMT Stationary Distribution (Waddington Landscape Basin)', 'FontSize', 16);

% Overlay the Bivalent / Hybrid boundaries from the Zhao et al. paper
hold on;
x_line = 0:0.1:N;
plot(x_line, x_line/2, 'w--', 'LineWidth', 2); % Upper bound of Activating
plot(x_line, 2*x_line, 'w--', 'LineWidth', 2); % Lower bound of Repressive
theta = linspace(0, pi/2, 100);
plot((N/4)*cos(theta), (N/4)*sin(theta), 'w--', 'LineWidth', 2); % Unmodified threshold
text(N/2, N/2, 'Hybrid E/M Zone', 'Color', 'white', 'FontSize', 12, 'FontWeight', 'bold', 'HorizontalAlignment', 'center');
hold off;

disp('Stationary distribution plotted.');