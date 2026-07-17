% run_master_bottleneck.m
clear all;

%% 1. Parameters (The Clinical Edge State)
N = 60; 
MAX_ITER = 6e4; 
period = 20; 
NS = 7; local = 1; division = 1; NSA = (NS-1)/2; 

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
t_drug_start = 250;   % HARDCODE: Drug hits exactly at t = 250
t_drug_end = 350;     % HARDCODE: Chemotherapy turns OFF (100 cell cycle pulse)
death_times = inf(num_runs, 1);
max_times = zeros(num_runs, 1);
all_R_trajectories = cell(num_runs, 1);
all_t_grids = cell(num_runs, 1);

disp('Simulating unified tumor ensemble and applying chemotherapy filter...');

%% 3. Simulate and Record (Executes only ONCE)
for i = 1:num_runs
    X0 = 1 * ones(N, 1); % Start 100% Epithelial
    
    [t_grid, A_interp, R_interp] = main_SSA_ensemble(@fnl, ka, kr, da, dr, eka, ekr, eda, edr, rho, ...
                                                     X0, MAX_ITER, period, N, NS, local, division);
    
    current_steps = length(t_grid);
    max_times(i) = t_grid(end);
    all_R_trajectories{i} = R_interp;
    all_t_grids{i} = t_grid;
    
       
    % Find exact drug start for this run
    this_drug_idx = find(t_grid >= t_drug_start, 1);
    
    if ~isempty(this_drug_idx)
        for j = this_drug_idx:current_steps
            if t_grid(j) > t_drug_end
                break; 
            end

            if R_interp(j) > death_threshold
                death_times(i) = t_grid(j); % Record exact death time
                break;
            end
        end
    end
    
    if mod(i, 10) == 0
        disp(['Completed ', num2str(i), '/', num2str(num_runs), ' cells.']);
    end
end

disp('Simulation complete. Generating linked figures...');

%% 4. FIGURE 1: The Population Survival Curve
plot_time_end = min(max_times);
master_t = linspace(0, plot_time_end, 1000);
survival_curve = zeros(1, length(master_t));

for k = 1:length(master_t)
    survival_curve(k) = sum(death_times > master_t(k));
end

figure('Name', 'Fig 1: Population Survival', 'Position', [100, 500, 700, 400], 'Color', 'w');
plot(master_t, survival_curve, 'k-', 'LineWidth', 2.5);
hold on;
xline(t_drug_start, 'r--', 'Chemotherapy Applied', 'LineWidth', 1.5, ...
      'LabelVerticalAlignment', 'bottom', 'FontSize', 12, 'FontWeight', 'bold');
xline(t_drug_end, 'b--', 'Chemotherapy Ended', 'LineWidth', 1.5, ...
    'LabelVerticalAlignment', 'top', 'FontSize', 12, 'FontWeight', 'bold');

xlim([0 plot_time_end]);
ylim([0 num_runs + 5]);
xlabel('Time (in Cell Cycles)', 'FontSize', 14, 'FontWeight', 'bold');
ylabel('Number of Surviving Cells', 'FontSize', 14, 'FontWeight', 'bold');
title('Tumor Population Dynamics and Persister Survival', 'FontSize', 16,'Color','k');
set(gca, 'XColor', 'k', 'YColor', 'k'); % Changes axis ticks/labels to black
cb.Color = 'k';                         % Changes colorbar ticks/border to black
cb.Label.Color = 'k';                   % Changes colorbar text label to black

grid on; box on; hold off;

%% 5. FIGURE 2: The Retroactive Trajectory Proof
figure('Name', 'Fig 2: Trajectory Proof', 'Position', [850, 500, 800, 400], 'Color', 'w');
hold on;

% Light background zones
yregion(40, 60, 'FaceColor', [1 0.5 0.5], 'FaceAlpha', 0.08); % Red Death Zone
yregion(0, 40, 'FaceColor', [0.5 0.7 1], 'FaceAlpha', 0.08);  % Blue Safe Zone

survivor_indices = find(death_times == inf);
doomed_indices = find(death_times ~= inf);

% BACKGROUND: Doomed cells (Faint Red)
for idx = 1:length(doomed_indices)
    i = doomed_indices(idx);
    t = all_t_grids{i};
    R = all_R_trajectories{i};
    death_idx = find(t >= death_times(i), 1);
    plot(t(1:death_idx), R(1:death_idx), 'Color', [0.8500 0.3250 0.0980 0.15], 'LineWidth', 0.5);
end

% FOREGROUND: Surviving cells (Bold Blue)
for idx = 1:length(survivor_indices)
    i = survivor_indices(idx);
    t = all_t_grids{i};
    R = all_R_trajectories{i};
    plot(t, R, 'Color', [0 0.4470 0.7410 0.9], 'LineWidth', 2.0);
end

xline(t_drug_start, 'k--', 'Chemotherapy Applied', 'LineWidth', 2.5, ...
      'LabelVerticalAlignment', 'top', 'FontSize', 14, 'FontWeight', 'bold');
xline(t_drug_end, 'k--', 'Chemotherapy Ended', 'LineWidth', 1.5, ...
    'LabelVerticalAlignment', 'top', 'FontSize', 14, 'FontWeight', 'bold');
yline(death_threshold, 'r-', 'Death Threshold', 'LineWidth', 2, 'LabelHorizontalAlignment', 'left');

xlim([0 min(cellfun(@max, all_t_grids))]);
ylim([0 60]);
xlabel('Time (in Cell Cycles)', 'FontSize', 14, 'FontWeight', 'bold');
ylabel('H3K27me3 Levels (Repressive Marks)', 'FontSize', 14, 'FontWeight', 'bold');
title('Trajectory of tumour cells', 'FontSize', 16,'Color','k');

set(gca, 'XColor', 'k', 'YColor', 'k'); % Changes axis ticks/labels to black
cb.Color = 'k';                         % Changes colorbar ticks/border to black
cb.Label.Color = 'k';                   % Changes colorbar text label to black

grid on; box on; hold off;

%% 6. Categorize the Survivors (Console Output)
hybrid_count = 0;
mesenchymal_count = 0;
relapsed_count = 0;

for idx = 1:length(survivor_indices)
    i = survivor_indices(idx);
    final_R = all_R_trajectories{i}(end);
    
    if final_R <= 40 && final_R >= 10
        hybrid_count = hybrid_count + 1;
    elseif final_R < 10
        mesenchymal_count = mesenchymal_count + 1;
    elseif final_R > 40
        relapsed_count = relapsed_count + 1; % Cells that went back to Epithelial
    end
end

disp('---Breakdown at the end of the sim---');
disp(['Total Initial Cells: ', num2str(num_runs)]);
disp(['Cells Killed by Drug: ', num2str(length(doomed_indices))]);
disp(['Total Surviving Persisters: ', num2str(length(survivor_indices))]);
disp(['  -> Hybrid E/M : ', num2str(hybrid_count)]);
disp(['  -> Pure Mesenchymal: ', num2str(mesenchymal_count)]);
disp(['  -> Relapsed back to Epithelial: ', num2str(relapsed_count)]);