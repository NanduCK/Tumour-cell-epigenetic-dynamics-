% analyze_kinetics.m
clear all;
load('EMT_ensemble_data.mat');

num_runs = length(A_ensemble);
dt = T_ensemble{1}(2) - T_ensemble{1}(1); % Time step size

%% 1. Phenotypic Autocorrelation (Memory Half-Life)
max_lags = round(500 / dt); % Look up to 500 cell cycles ahead
C_tau_sum = zeros(2*max_lags + 1, 1);

for i = 1:num_runs
    % Force into column vectors
    A = A_ensemble{i}(:);
    R = R_ensemble{i}(:);
    
    % Convert to 1D Phenotypic axis (Active - Repressive)
    M = A - R;
    
    % Normalize by subtracting mean to get fluctuations
    M_norm = M - mean(M);
    
    % Compute autocorrelation
    [C, lags] = xcorr(M_norm, max_lags, 'normalized');
    
    % Ensure C is added as a column vector
    C_tau_sum = C_tau_sum + C(:);
end

C_tau_avg = C_tau_sum / num_runs;
tau_lags = lags * dt;

% Plot Autocorrelation
figure('Position', [100, 100, 1000, 400]);
subplot(1,2,1);
plot(tau_lags(tau_lags >= 0), C_tau_avg(tau_lags >= 0), 'k-', 'LineWidth', 2);
xlabel('Time Lag \tau (Cell Cycles)', 'FontSize', 12);
ylabel('Autocorrelation C(\tau)', 'FontSize', 12);
title('Epigenetic Memory Decay', 'FontSize', 14);
grid on;
xlim([0 200]); % Zoom in on the decay

%% 2. Hybrid E/M Residence Time (Mean Dwell Time)
dwell_times = [];

for i = 1:num_runs
    % Force into column vectors
    A = A_ensemble{i}(:);
    R = R_ensemble{i}(:);
    T = T_ensemble{i}(:);
    
    % Define the Bivalent (Hybrid E/M) mathematical conditions:
    is_hybrid = (A.^2 + R.^2 > (N/4)^2) & (A <= 2.*R) & (R <= 2.*A);
    
    % Find entry and exit points (now strictly column vectors)
    transitions = diff([0; is_hybrid; 0]); 
    entry_indices = find(transitions == 1);
    exit_indices = find(transitions == -1) - 1;
    
    for j = 1:length(entry_indices)
        duration = T(exit_indices(j)) - T(entry_indices(j));
        % Ignore micro-fluctuations (e.g., staying in state for less than 1 cell cycle)
        if duration > 1 
            dwell_times = [dwell_times; duration];
        end
    end
end

% Plot Residence Time Distribution
subplot(1,2,2);
histogram(dwell_times, 50, 'Normalization', 'pdf', 'FaceColor', [0.8 0.2 0.2]);
xlabel('Dwell Time (Cell Cycles)', 'FontSize', 12);
ylabel('Probability Density', 'FontSize', 12);
mean_dwell = mean(dwell_times);
title(['Hybrid E/M Residence Time (Mean: ', num2str(round(mean_dwell,1)), ')'], 'FontSize', 14);
grid on;

disp(['Mean Hybrid E/M Dwell Time: ', num2str(mean_dwell), ' cell cycles.']);