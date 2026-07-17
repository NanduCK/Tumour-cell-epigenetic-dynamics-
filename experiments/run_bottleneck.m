% run_bottleneck_standalone.m
clear all;

%% 1. Parameters (The Symmetric Bivalent State)
N = 60; 
MAX_ITER = 3e4; 
period = 20; 
NS = 7; local = 0; division = 1; NSA = (NS-1)/2; 

gamma = 0.2;
ka1 = 2.0; % Activating drive
kr1 = 2.5; % Symmetric Polycomb drive
da = 1 * ones(NSA,1); eda = gamma * da;   
dr = 1 * ones(NSA,1); edr = gamma * dr;  
ka = ka1 * ones(NSA,1); eka = gamma * ka; 
kr = kr1 * ones(NSA,1); ekr = gamma * kr;
rho = 0; 

num_runs = 100;

%% 2. Drug Bottleneck Setup
death_threshold = 40; % >40 H3K27me3 marks = Epithelial = Death
t_drug_start = 0; % Will be set dynamically during Run 1
death_times = inf(num_runs, 1); % Record when each cell dies (Inf = survived)
max_times = zeros(num_runs, 1); % Record the physical end time of each simulation

disp('Simulating tumor ensemble and applying chemotherapy filter...');
tic;

%% 3. Simulate and Filter Simultaneously
for i = 1:num_runs
    X0 = 1 * ones(N, 1); 
    
    % Generate the trajectory for Cell i
    [t_grid, A_interp, R_interp] = main_SSA_ensemble(@fnl, ka, kr, da, dr, eka, ekr, eda, edr, rho, ...
                                                     X0, MAX_ITER, period, N, NS, local, division);
    
    current_steps = length(t_grid);
    max_times(i) = t_grid(end);
    
    % Set drug application time based on the first run (at 20% of its physical time)
    if i == 1
        t_drug_start = t_grid(round(current_steps * 0.5));
        disp(['Drug will be applied at physical time t = ', num2str(t_drug_start)]);
    end
    
    % Find the array index where the drug phase starts for THIS specific run
    drug_start_index = find(t_grid >= t_drug_start, 1);
    
    % Check for death condition strictly after drug administration
    if ~isempty(drug_start_index)
        for j = drug_start_index:current_steps
            if R_interp(j) > death_threshold
                % The cell fell into the sensitive Epithelial basin and died.
                death_times(i) = t_grid(j); % Record the exact physical time of death
                break; % Move to the next cell in the ensemble
            end
        end
    end
    
    % Progress update
    if mod(i, 10) == 0
        disp(['Completed ', num2str(i), '/', num2str(num_runs), ' cells.']);
    end
end

toc;

%% 4. Process Survival Data and Plot
% Create a universal time grid for plotting, ending at the shortest simulation time
% (to ensure we don't artificially report cells dying just because their simulation ended early)
plot_time_end = min(max_times);
master_t = linspace(0, plot_time_end, 1000);
survival_curve = zeros(1, length(master_t));

% Calculate how many cells are alive at each point in time
for k = 1:length(master_t)
    current_time = master_t(k);
    % A cell is alive if its death time is strictly greater than the current time
    survival_curve(k) = sum(death_times > current_time);
end

surviving_cells = sum(death_times > plot_time_end);
disp(['Final Surviving Persister Cells: ', num2str(surviving_cells)]);

% Plotting
figure('Position', [100, 100, 700, 400], 'Color', 'w');
plot(master_t, survival_curve, 'k-', 'LineWidth', 2.5);
hold on;

% Draw the drug administration line
xline(t_drug_start, 'r--', 'Chemotherapy Applied', 'LineWidth', 1.5, ...
      'LabelVerticalAlignment', 'bottom', 'FontSize', 12, 'FontWeight', 'bold');

% Formatting
xlim([0 plot_time_end]);
ylim([0 num_runs + 5]);
xlabel('Time (in Cell Cycles)', 'FontSize', 14, 'FontWeight', 'bold');
ylabel('Number of Surviving Cells', 'FontSize', 14, 'FontWeight', 'bold');
title('Tumor Population Collapse and Persister Survival', 'FontSize', 16);
grid on; box on; hold off;