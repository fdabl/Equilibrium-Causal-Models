# ################################################################
# This script analyzes the output of simulation-backshift-create.R
# ################################################################
source('ECM.R')
source('helpers.R')

cores <- detectCores()
cl <- makeCluster(cores - 8)
registerDoParallel(cl)

B4 <- t(matrix(
  c(.5,   0, 0,  0,
    .3, .4, .2,  0,
    0,  .2, .2, -.4,
    0,   0,  0,  .4), 4, 4, byrow = T
))

Bs <- list('4' = t(rs_beta(t(B4))))
res <- readRDS('Results/backshift_estimates.RDS')

start <- Sys.time()

# Analyze the estimated graphs
dat <- foreach(i = seq(nrow(res)), .combine = 'rbind') %dopar% {
  
  cur <- res[i, ]
  Bhat <- cur$estimate$Ahat
  p <- ncol(Bhat)
  Btrue <- Bs[[as.character(p)]]
  Sigma <- cur$Sigmahat
  
  diag(Bhat) <- NA
  diag(Btrue) <- NA
  
  metrics <- get_met(Bhat, Btrue, Sigma, Sigma)
  
  c(cur$config, metrics)
}

colnames(dat) <- c(
  'times', 'n', 'p', 'corr', 'nr_int', 'nr_targets',
  'effect_size', 'warning', 'L1', 'L2', 'mult_median',
  'mult_mean', 'mult_min', 'mult_max',
  paste0('L1_X', seq(4)),
  paste0('Mult_median_X', seq(4)),
  paste0('Mult_mean_X', seq(4)),
  paste0('Mult_min_X', seq(4)),
  paste0('Mult_max_X', seq(4))
)

dat <- data.frame(dat)
dat <- apply(dat, 2, as.numeric)
rownames(dat) <- NULL
write.csv(dat, 'Results/backshift_metrics_test.csv', row.names = FALSE)

end <- Sys.time()
print(end - start)
