# Tumour-cell-epigenetic-dynamics

#Stochastic Modelling of EMT and Drug Resistance

This repository contains a MATLAB-based computational framework for simulating the dynamics of bivalent epigenetic domains at the Epithelial-Mesenchymal Transition (EMT) loci in tumour cells. 

Building upon the core simulation engine adapted from Zhao et al. (2021), this project was done as part of the PHCCO program at IISc, Bangalore and utilises the Gillespie Stochastic Simulation Algorithm (SSA) to model chromatin state transitions (H3K4me3 vs H3K27me3). The analysis pipeline specifically investigates epigenetic memory, Waddington landscape stationary distributions, and how transient Hybrid E/M states act as survival bottlenecks during acute chemotherapy exposure.

## 🗂️ Repository Structure

The codebase is modularised into three distinct functions:

```text
├── src/                        # Core mathematical engine
│   ├── fl.m                    # Linear Feedback condition
│   ├── fnl.m                   # Non-linear Feedback     
│   └── main_SSA_ensemble.m     # Gillespie (SSA) engine for chromatin transitions
├── experiments/                # Controllers for generating simulation data
│   ├── run_ensemble.m          # Base ensemble generation
│   ├── run_sweep.m             # Parameter sweeping for repressive methylation
│   ├── run_master_dtp.m        # Population dynamics simulations
│   └── run_trajectory_proof.m  # Visual tracing of persister cell origins
└── analysis/                   # Analytics and visualisation scripts
    ├── analyze_kinetics.m      # Autocorrelation and residence time mapping
    ├── analyze_stationary.m    # 2D Waddington landscape density plotting
    └── plot_bifurcation.m      # Stochastic bifurcation diagram generation
```
