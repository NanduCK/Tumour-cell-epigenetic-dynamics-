% run_trajectory_proof.m
clear all;

%% 1. Parameters (The Clinical Edge State)
N = 60; 
MAX_ITER = 3e4; 
period = 20; 
NS = 7; local = 0; division = 1; NSA = (NS-1)/2; 

gamma = 0.2;
ka1 = 2.0; 
kr1 = 2.5; % The deep Epithelial basin
da = 1 * ones(NSA,1); eda = gamma * da;   
dr = 1 * ones(NSA,1); edr = gamma * dr;  
ka = ka1 * ones(NSA,1); eka = gamma * ka; 
kr = kr1 * ones(NSA,1); ekr = gamma * kr;
rho = 0; 

num_runs = 100;

%% 2. Drug Bottleneck Setup
death_threshold = 40; % >40 H3K27me3 marks = Epithelial = Death
death_times = inf(num_runs, 1);
all_R_trajectories = cell(num_runs, 1);
all_t_grids = cell(num_runs, 1);

disp('Simulating trajectories and sorting by survival fate...');

%% 3. Simulate and Record Every Trajectory
for i = 1:num_runs
    X0 = 1 * ones(N, 1); % Start 100% Epithelial

    [t_grid, A_interp, R_interp] = main_SSA_ensemble(@fnl, ka, kr, da, dr, eka, ekr, eda, edr, rho, ...
        X0, MAX_ITER, period, N, NS, local, division);

    current_steps = length(t_grid);
    all_R_trajectories{i} = R_interp;
    all_t_grids{i} = t_grid;

    if i == 1
        drug_start_index = round(current_steps * 0.5);
        t_drug_start = t_grid(drug_start_index);
    end

    % Find exact drug start for this run
    this_drug_idx = find(t_grid >= t_drug_start, 1);

    if ~isempty(this_drug_idx)
        for j = this_drug_idx:current_steps
            if R_interp(j) > death_threshold
                death_times(i) = t_grid(j); % Record exact death time
                break;
            end
        end
    end
end

disp(['Simulation complete. Plotting visual proof...']);

%% 4. Plot the Trajectory Proof (Clean Layered Rendering)
figure('Position', [100, 100, 800, 500], 'Color', 'w');
hold on;

% 1. Lighten the background zones so they don't fight the data
yregion(40, 60, 'FaceColor', [1 0.5 0.5], 'FaceAlpha', 0.08); % Ultra-faint red
yregion(0, 40, 'FaceColor', [0.5 0.7 1], 'FaceAlpha', 0.08);  % Ultra-faint blue

% 2. Sort the cells by their fate (Z-Ordering preparation)
survivor_indices = find(death_times == inf);
doomed_indices = find(death_times ~= inf);

% 3. BACKGROUND LAYER: Plot DOOMED cells first
for idx = 1:length(doomed_indices)
    i = doomed_indices(idx);
    t = all_t_grids{i};
    R = all_R_trajectories{i};
    death_idx = find(t >= death_times(i), 1);
    
    % Use RGBA color: [Red Green Blue Alpha]. Alpha = 0.15 makes it mostly transparent.
    plot(t(1:death_idx), R(1:death_idx), 'Color', [0.8500 0.3250 0.0980 0.15], 'LineWidth', 0.5);
end

% 4. FOREGROUND LAYER: Plot SURVIVING cells last
for idx = 1:length(survivor_indices)
    i = survivor_indices(idx);
    t = all_t_grids{i};
    R = all_R_trajectories{i};
    
    % Use bold, fully opaque blue to pop to the front
    plot(t, R, 'Color', [0 0.4470 0.7410 0.9], 'LineWidth', 2.0);
end

% 5. Add Annotations
xline(t_drug_start, 'k--', 'Chemotherapy Applied', 'LineWidth', 2.5, ...
      'LabelVerticalAlignment', 'top', 'FontSize', 14, 'FontWeight', 'bold');

yline(death_threshold, 'r-', 'Death Threshold', 'LineWidth', 2, 'LabelHorizontalAlignment', 'left');

% Formatting
xlim([0 min(cellfun(@max, all_t_grids))]);
ylim([0 60]);
xlabel('Time (in Cell Cycles)', 'FontSize', 14, 'FontWeight', 'bold');
ylabel('H3K27me3 Levels (Repressive Marks)', 'FontSize', 14, 'FontWeight', 'bold');
title('Retroactive Trajectory Proof of Persister Origins', 'FontSize', 16, Color='k');

set(gca, 'XColor', 'k', 'YColor', 'k'); % Changes axis ticks/labels to black
cb.Color = 'k';                         % Changes colorbar ticks/border to black
cb.Label.Color = 'k';                   % Changes colorbar text label to black

grid on; box on; hold off;