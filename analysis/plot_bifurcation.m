% plot_bifurcation.m
clear all;
load('EMT_sweep_data.mat');

% Prepare grid for plotting
% X-axis: kr values (Parameter)
% Y-axis: Net Phenotype (-N to +N)
X = kr_values;
Y = phenotype_bins(1:end-1) + 0.5; % Bin centers
Z = sweep_density'; % Transpose so Y is Phenotype and X is Parameter

%% Plot the Stochastic Bifurcation
figure('Position', [100, 100, 800, 600], 'Color', 'w');

% Use imagesc or pcolor to create the density heatmap
imagesc(X, Y, Z);
set(gca, 'YDir', 'normal'); % Keep Mesenchymal (+N) at the top, Epithelial (-N) at the bottom
colormap(parula); % Parula or 'jet' work well for this

% Formatting
cb = colorbar;
cb.Label.String = 'Probability Density';
cb.Label.FontSize = 12;

xlabel('Polycomb Recruitment Rate (kr)', 'FontSize', 14, 'FontWeight', 'bold');
ylabel('Net Phenotype (H3K4me3 - H3K27me3)', 'FontSize', 14, 'FontWeight', 'bold');
title('Stochastic Bifurcation Diagram: EMT Plasticity', 'FontSize', 16, 'Color','k');

set(gca, 'XColor', 'k', 'YColor', 'k'); % Changes axis ticks/labels to black
cb.Color = 'k';                         % Changes colorbar ticks/border to black
cb.Label.Color = 'k';                   % Changes colorbar text label to black

% Overlay guidelines
hold on;
yline(0, 'w--', 'LineWidth', 1.5); % Center Hybrid E/M line
xline(2.0, 'r-', 'LineWidth', 1.5); % The symmetric threshold (ka = 2.0)
text(2.05, 50, '\leftarrow Symmetric Point (ka = kr)', 'Color', 'red', 'FontSize', 12);
hold off;

% Adjust limits
ylim([-N N]);
xlim([min(kr_values) max(kr_values)]);

% Optional: Set log scale if you want to highlight the shallow basins 
% set(gca, 'ColorScale', 'log');