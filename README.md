# Target-trial-emulation
This folder contains R code used to produce the results for the paper:
“Liberal or restrictive transfusion for veno-arterial Extracorporeal Membrane Oxygenation patients: a target trial emulation using the OBLEX study data” — Thao Le et al.

There are three R scripts included:
- 01_AuxFun.R: Contains auxiliary functions used throughout the analysis.
- 02_Analysis.R: Runs the main analysis for the paper:
	- Generates sequential trials at each time point
	- Calculates weights
	- Fits an Aalen additive hazards model
	- Estimates survival probabilities for each transfusion strategy
	- Generate figure 2B: Estimated survival difference between the liberal and restrictive transfusion strategies with bootstrap confidence intervals
- 03_BootstrapCI.R: Computes confidence intervals for the estimated survival probabilities using bootstrap sampling.
