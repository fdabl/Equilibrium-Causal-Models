source('ECM.R')
source('helpers.R')


B4 <- t(matrix(
  c(.5,   0, 0,  0,
    .3, .4, .2,  0,
    0,  .2, .2, -.4,
    0,   0,  0,  .4), 4, 4, byrow = T
))

B <- rs_beta(t(B4))


########################################
# The following results are for Figure 7
########################################
p <- 4
ns <- 500
times <- 250
true_trait_prop <- 0.70

trait_prop <- seq(0.50, 0.90, 0.025)
configs <- expand.grid(
  n = ns, time = seq(times),
  trait_prop = trait_prop
)
nr_config <- nrow(configs)

set.seed(1)
registerDoParallel(cores = 9)

# Estimate B from equilibrium and snapshot data
res_sem <- foreach(i = seq(nr_config), .combine = 'rbind') %dopar% {
  
  config <- configs[i, ]
  n <- config$n
  trait_prop <- config$trait_prop
  
  # Create data with certain trait / state variance composition
  dat <- get_measurement_dat(n = n, trait_prop = true_trait_prop)
  
  matrices <- get_lavaan_measurement(dat, trait_prop = trait_prop)
  Bhat <- t(matrices$Bhat)
  Sigmahat <- matrices$Sigmahat
  inthat <- rep(0, p)
  
  press_X1 <- comp_press(t(Bhat), inthat, target = 1, a = 1)
  press_X4 <- comp_press(t(Bhat), inthat, target = 4, a = 1)
  beta_names <- c('B12', 'B23', 'B32', 'B43')
  beta_values <- c(Bhat[1, 2], Bhat[2, 3], Bhat[3, 2], Bhat[4, 3])
  
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
  dat$trait_prop <- trait_prop
  dat$true_trait_prop <- true_trait_prop
  
  list(
    'dat' = dat,
    'Bhat' = Bhat,
    'Sigma_hat' = Sigmahat
  )
}

saveRDS(res_sem, 'Results/state_trait_results.RDS')


##########################################
# The following results are for Appendix E
##########################################

### Simulate data
ns <- 2000
times <- 250

resid_prop <- seq(0, 0.50, 0.10)
resid_corrs <- c(-0.25, 0, 0.25)
configs <- expand.grid(
  n = ns, time = seq(times),
  resid_prop = resid_prop, resid_cor = resid_corrs
)
nr_config <- nrow(configs)

set.seed(1)
registerDoParallel(cores = 9)

res_sem <- foreach(i = seq(nr_config), .combine = 'rbind') %dopar% {
  
  config <- configs[i, ]
  n <- config$n
  resid_cor <- config$resid_cor
  resid_prop <- config$resid_prop
  
  # Create data with certain trait / state variance composition
  dat <- get_measurement_dat(n = n, trait_prop = 1 - resid_prop, resid_cor = resid_cor)
  matrices <- get_lavaan_measurement(dat, trait_prop = 1 - resid_prop)
  Bhat <- t(matrices$Bhat)
  Sigmahat <- matrices$Sigmahat
  inthat <- rep(0, p)
  
  cnames <- paste0('X', seq(p))
  press_X1 <- comp_press(t(Bhat), inthat, target = 1, a = 1)
  press_X4 <- comp_press(t(Bhat), inthat, target = 4, a = 1)
  rownames(press_X1) <- rownames(press_X4) <- cnames
  
  dat <- data.frame(
    rbind(
      t(press_X1),
      t(press_X4)
    )
  )
  
  dat$target <- factor(ifelse(dat$X1 == 1, 'X1', 'X4'))
  dat <- reshape2::melt(dat, c('X1', 'X4', 'target'))
  dat$time <- config$time
  dat$resid_prop <- resid_prop
  dat$resid_cor <- resid_cor
  
  list(
    'dat' = dat,
    'Bhat' = Bhat,
    'Sigma_hat' = Sigmahat,
    'resid_cor' = resid_cor,
    'resid_prop' = resid_prop
  )
}

saveRDS(res_sem, 'Results/measurement_corr_results.RDS')
