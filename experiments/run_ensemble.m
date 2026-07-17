% run_ensemble.m
clear all;

%% 1. Simulation Parameters 
num_runs = 100;      
MAX_ITER = 2e5;     
period = 20; 
N = 60;              
NS = 7;
local = 0;          
division = 1; 
NSA = (NS-1)/2; 

ka1 = 3; kr1 = 3; gamma = 0.2;
da = 1 * ones(NSA,1); eda = gamma * da;   
dr = 1 * ones(NSA,1); edr = gamma * dr;  
ka = ka1 * ones(NSA,1); eka = gamma * ka; 
kr = kr1 * ones(NSA,1); ekr = gamma * kr;
rho = 0; 

%% 2. Preallocate Cell Arrays
A_ensemble = cell(num_runs, 1);
R_ensemble = cell(num_runs, 1);
T_ensemble = cell(num_runs, 1);

disp(['Starting sequential ensemble simulation of ', num2str(num_runs), ' runs...']);
tic;

%% 3. The Sequential Loop
for i = 1:num_runs
    % X0 initialized INSIDE the loop so each run gets a fresh state
    X0 = 3 * ones(N, 1); 
    
    % Run the modified Gillespie simulation
    [t_grid, A_interp, R_interp] = main_SSA_ensemble(@fnl, ka, kr, da, dr, eka, ekr, eda, edr, rho, ...
                                                     X0, MAX_ITER, period, N, NS, local, division);
    
    % Store the standardized results
    T_ensemble{i} = t_grid;
    A_ensemble{i} = A_interp;
    R_ensemble{i} = R_interp;
    
    % Optional: Display progress so you know it hasn't frozen
    if mod(i, 2) == 0
        disp(['Completed run ', num2str(i), ' of ', num2str(num_runs)]);
    end
end

execution_time = toc;
disp(['Ensemble completed in ', num2str(execution_time), ' seconds.']);

%% 4. Save the Data
save('EMT_ensemble_data.mat', 'T_ensemble', 'A_ensemble', 'R_ensemble', 'N');
disp('Data saved to EMT_ensemble_data.mat');