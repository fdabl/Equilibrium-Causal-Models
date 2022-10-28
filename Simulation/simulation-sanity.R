# ###########################################################################
# This script runs the "sanity check" simulation study summarized in Figure 4
# ###########################################################################
source('functions_helpers.R')
source('functions_ecm.R')


# True underlying p = 4 system
B4 <- t(matrix(
  c(.5,   0, 0,  0,
    .3, .4, .2,  0,
    0,  .2, .2, -.4,
    0,   0,  0,  .4), 4, 4, byrow = T
))

# Equilibrium matrix
B <- t(rs_beta(t(B4)))

### Simulate data
p <- ncol(B)
ns <- seq(50, 1000, 50)
times <- 250
configs <- expand.grid(n = ns, time = seq(times))
nr_config <- nrow(configs)

set.seed(1)
registerDoParallel(cores = 9)

# Estimate B from equilibrium and snapshot data
res_sem <- foreach(i = seq(nr_config), .combine = 'rbind') %dopar% {
  
  config <- configs[i, ]
  n <- config$n
  
  # Draw intercepts and create equilibrium data
  intercepts <- rmvnorm(n, mean = rep(0, p), sigma = diag(p))
  dat_eq <- create_eqdat(B, intercepts)
  
  # Snapshot data
  # Note that estimates get very large at some point with high noise variance
  matrices_eq <- get_lavaan(dat_eq, t(B), estimate_sigma = FALSE)
  Beq <- t(matrices_eq$Bhat)
  Sigmaeq <- matrices_eq$Sigmahat
  inteq <- matrices_eq$intercepthat
  
  press_X1 <- comp_press(t(Beq), inteq, target = 1, a = 1)
  press_X4 <- comp_press(t(Beq), inteq, target = 4, a = 1)
  beta_names <- c('B12', 'B23', 'B32', 'B43')
  beta_values <- c(Beq[1, 2], Beq[2, 3], Beq[3, 2], Beq[4, 3])
  
  causal_names <- c('C12', 'C13', 'C43', 'C42')
  causal_values <- c(press_X1[2], press_X1[3], press_X4[3], press_X4[2])
  
  dat <- data.frame(
    betas = beta_names,
    causal = causal_names,
    beta_values = beta_values,
    causal_values = causal_values
  )
  
  dat$time <- config$time
  dat$n <- n
  
  list(
    'dat' = dat,
    'Bhat' = Beq,
    'Sigma_hat' = Sigmaeq
  )
}

saveRDS(res_sem, 'Results/sanity-results.RDS')
