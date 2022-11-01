# #########################################################################
# This script runs the measurement simulation study summarized in Figure 6
# #########################################################################
source('ECM.R')
source('helpers.R')


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
ns <- 2000
times <- 250

noise_prop <- seq(0, 0.50, 0.10)
resid_vars <- noise_prop / (1 - noise_prop)
resid_corrs <- c(-0.25, 0, 0.25)
configs <- expand.grid(
  n = ns, time = seq(times),
  resid_var = resid_vars, resid_cor = resid_corrs
)
nr_config <- nrow(configs)

set.seed(1)
registerDoParallel(cores = 9)

# Estimate B from equilibrium and snapshot data
res_sem <- foreach(i = seq(nr_config), .combine = 'rbind') %dopar% {
  
  config <- configs[i, ]
  n <- config$n
  resid_var <- config$resid_var
  resid_cor <- config$resid_cor
  
  Sigma <- diag(resid_var, nrow = p, ncol = p)
  Sigma[1, 2] <- Sigma[2, 1] <- (resid_cor * resid_var)
  
  # Draw intercepts and create equilibrium data
  intercepts <- rmvnorm(n, mean = rep(0, p), sigma = diag(p))
  dat_eq <- create_eqdat(B, intercepts)
  dat_noise <- create_noisedat(B, intercepts, Sigma)
  dat_snap <- dat_eq + dat_noise
  
  # Snapshot data
  # Note that estimates get very large at some point with high noise variance
  matrices_snap <- get_lavaan(dat_snap, t(B), estimate_sigma = FALSE)
  Bsnap <- t(matrices_snap$Bhat)
  Sigmasnap <- matrices_snap$Sigmahat
  intsnap <- matrices_snap$intercepthat
  
  press_snap_X1 <- comp_press(t(Bsnap), intsnap, target = 1, a = 1)
  press_snap_X4 <- comp_press(t(Bsnap), intsnap, target = 4, a = 1)
  
  dat_int <- data.frame(
    rbind(
      t(press_snap_X1),
      t(press_snap_X4)
    )
  )
  resid_prop <- resid_var / (1 + resid_var)
  dat_int$target <- factor(ifelse(dat_int$X1 == 1, 'X1', 'X4'))
  dat_int <- reshape2::melt(dat_int, c('X1', 'X4', 'target'))
  dat_int$time <- config$time
  dat_int$resid_var <- resid_var
  dat_int$resid_prop <- resid_prop
  dat_int$resid_cor <- resid_cor
  
  list(
    'dat' = dat_int,
    'Bhat' = Bsnap,
    'Sigma_hat' = Sigmasnap,
    'resid_cor' = resid_cor,
    'resid_var' = resid_var,
    'resid_prop' = resid_prop
  )
}

saveRDS(res_sem, 'Results/measurement_results.RDS')
