# Example to show how to estimate an equilibrium causal model using Backshift
# Data are assumed to consist of equilibrium measurements
source('helpers.R')

p <- 4

# Specify beta matrix
a <- 0.5
b <- 0.33
c <- 0.25
d <- -0.5

B <- matrix(c(0,0,0,0,
              a,0,b,0,
              0,c,0,d,
              0,0,0,0),4,4,byrow = T)

n <- 1000
vars <- seq(p)
effect_size <- 1
nr_targets <- 2
nr_interventions <- 3

# For each intervention, randomly select the variables to intervene on
ints <- list()
for (j in seq(nr_interventions)) {
  target <- sample(vars, 3)
  int <- rep(0, p)
  int[target] <- 1
  
  ints[[j]] <- int
}

ints[[j + 1]] <- rep(0, p) # observational setting

# Create error covariance and data
Sigma <- create_Sigma(p, noise_sd = 1, corr = 0)
backshift_dat <- create_interventions_data(
  B, Sigma, n, ints, intervention_sd = 1, effect_size = effect_size
)

# indicator variable indicates which portion of the data
# belongs to which environment, needed for backShift
nr_environments <- length(ints)
indicator <- rep(seq(nr_environments), each = n)

# See ?backShift for details
res <- backShift(backshift_dat, indicator, threshold = 0.75, ev = 1, nsim = 500)

# Compare estimate to ground truth
B
res$Ahat
res$Ahat * ifelse(res$AhatAdjacency, 1, 0) # Only select edges that survive stability selection
