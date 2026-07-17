% run_sweep.m
clear all;

%% 1. Sweep Parameters
% Sweeping the repressive recruitment rate (Polycomb drive)
% We sweep from 1.0 (highly deficient/plastic) to 3.0 (highly rigid/epithelial)
kr_values = linspace(1.0, 3.0, 21); 
num_points = length(kr_values);

%% 2. Gillespie Engine Parameters
num_runs = 50;      % 50 runs per parameter point is usually enough for a clean density
MAX_ITER = 1e5;     % Keep it relatively short to save sweep time
period = 20; 
N = 60;             % System size 
NS = 7;
local = 0;          % Global interactions 
division = 1; 
NSA = (NS-1)/2; 

ka1 = 2.0;          % Activating drive remains locked at 2.0
gamma = 0.2;
da = 1 * ones(NSA,1); eda = gamma * da;   
dr_base = 1;        % Basal demethylation rate 

%% 3. Preallocate Sweep Data Storage
% We will store the 1D Net Phenotype (Active - Repressive) histogram for each point
phenotype_bins = -N:N; 
sweep_density = zeros(num_points, length(phenotype_bins)-1);

disp('Starting Stochastic Parameter Sweep...');
tic;

%% 4. The Sweep Loop
for p = 1:num_points
    kr1 = kr_values(p);

    % Update dynamic arrays for this sweep point
    dr = dr_base * ones(NSA,1); edr = gamma * dr;  
    ka = ka1 * ones(NSA,1); eka = gamma * ka; 
    kr = kr1 * ones(NSA,1); ekr = gamma * kr;
    rho = 0; 

    A_flat = [];
    R_flat = [];

    % Run the ensemble for the current parameter
    for i = 1:num_runs
        X0 = 3 * ones(N, 1); 
        [t_grid, A_interp, R_interp] = main_SSA_ensemble(@fnl, ka, kr, da, dr, eka, ekr, eda, edr, rho, ...
            X0, MAX_ITER, period, N, NS, local, division);

        % Discard burn-in (first 20%)
        start_idx = round(0.2 * length(t_grid)) + 1;
        A_flat = [A_flat; A_interp(start_idx:end)'];
        R_flat = [R_flat; R_interp(start_idx:end)'];
    end

    % Calculate Net Phenotype: H3K4me3 (Active) - H3K27me3 (Repressive)
    Net_Phenotype = A_flat - R_flat;

    % Bin the data to create the probability density for this parameter point
    [counts, ~] = histcounts(Net_Phenotype, phenotype_bins, 'Normalization', 'probability');
    sweep_density(p, :) = counts;

    disp(['Completed point ', num2str(p), '/', num2str(num_points), ' (kr = ', num2str(kr1), ')']);
end

execution_time = toc;
disp(['Sweep completed in ', num2str(execution_time/60), ' minutes.']);

save('EMT_sweep_data.mat', 'kr_values', 'phenotype_bins', 'sweep_density', 'N');