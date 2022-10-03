# Equilibrium Causal Models
This repository contains code to reproduce all analyses and figures in Ryan, O.<sup>&#11089;</sup> & Dablander, F.<sup>&#11089;</sup> ([2022](LINK)). Equilibrium Causal Models: Connecting Dynamical Systems Modeling and Cross-Sectional Data Analysis.

- **Code/** holds relevant R code:
    - **Code/helpers.R** holds useful R functions.
    - **Code/simulation-sanity.R** runs the simulation summarized in Figure 4.
    - **Code/simulation-measurement.R** runs the simulation summarized in Figure 6.
    - **Code/simulation-backshift-create.R** creates data and estimates causal models using backshift.
    - **Code/simulation-backshift-analyze.R** calculates metrics of estimated causal models.
    - **Code/figures.R** provides code to recreate all figures in the manuscript.
- **Results/** includes all simulation results discussed in the manuscript:
    - **Results/sanity-results.RDS** holds the results of the sanity check simulation.
    - **Results/measurement-results.RDS** holds the results of the measurement simulation.
    - **Results/backshift-estimates.RDS** holds the estimated causal models.
    - **Results/backshift-metrics.RDS** holds the metrics calculated on the estimated causal models.
- **Figures/** includes all figures in the manuscript.
