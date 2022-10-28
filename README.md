# Equilibrium Causal Models
This repository contains code to reproduce all analyses and figures in Ryan, O.<sup>&#11089;</sup> \& Dablander, F.<sup>&#11089;</sup> ([2022](LINK)). Equilibrium Causal Models: Connecting Dynamical Systems Modeling and Cross-Sectional Data Analysis.

- `figures.R` code to reproduce all figures in the manuscript
- `functions_ecm.R` helper functions relating to ECMs, such as computing shift and press interventions, re-scaling matrices, etc.
- `functions_helpers.R` useful R functions for simulation and plotting
- **Simulation/** holds R code to reproduce all simulation studies in the paper:
    - `Simulation/simulation-sanity.R` runs the simulation summarized in Figure 4.
    - `Simulation/simulation-measurement.R` runs the simulation summarized in Figure 6.
    - `Simulation/simulation-backshift-create.R` creates data and estimates causal models using backshift.
    - `Simulation/simulation-backshift-analyze.R` calculates metrics of estimated causal models (summarized in Figures 10 and 11)
- **Results/** includes all simulation results discussed in the manuscript:
    - `Results/sanity-results.RDS` holds the results of the sanity check simulation.
    - `Results/measurement-results.RDS` holds the results of the measurement simulation.
    - `Results/backshift-estimates.RDS` holds the estimated causal models.
    - `Results/backshift-metrics.RDS` holds the metrics calculated on the estimated causal models.
- **Example/** includes code for an example `lavaan` analysis based on simulated data
    - `Example/example_datagen.R` creates simulated data based on the running example in the main text with imperfect measurements of the equilibrium
    - `Example/example_data.RDS` simulated dataset created by `example_datagen.R`
    - `Example/example_modelfit.R` fits the ECM using lavaan and the measurement model constraints described in Appendix C. 
- **Figures/** includes all figures in the manuscript.
