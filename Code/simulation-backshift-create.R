# ########################################################################
# This script runs the backshift simulation study summarized in Appendix C
# ########################################################################
source('Code/helpers.R')

cores <- detectCores()
cl <- makeCluster(cores - 5)
registerDoParallel(cl)

B4 <- t(matrix(
  c(.5,   0, 0,  0,
    .3, .4, .2,  0,
    0,  .2, .2, -.4,
    0,   0,  0,  .4), 4, 4, byrow = T
))

ps <- 4
times <- seq(100)
corrs <- c(0, 0.50)
nr_targets <- seq(max(ps))
effect_sizes <- c(0.50, 1, 2)
nr_interventions <- seq(3, 10)
ns <- c(250, 500, 1000, 2500, 5000)
Bs <- list('4' = t(rs_beta(t(B4))))

configs <- expand.grid(
  times = times, n = ns, p = ps,
  nr_targets = nr_targets,
  nr_interventions = nr_interventions,
  corrs = corrs, effect_sizes = effect_sizes
)

configs <- configs[configs$p >= configs$nr_targets, ]
nr_config <- nrow(configs)

start <- Sys.time()

# Create data and estimate causal model using backshift
res <- foreach(
  i = seq(nr_config), .combine = 'rbind',
  .export = c(
    'backShift', 'BACKSHIFT', 'create_interventions_data',
    'create_Sigma', 'is.positive.definite', 'rmvnorm'
    )
  ) %dopar% {

  conf <- configs[i, ]
  time <- conf$times
  p <- conf$p
  n <- conf$n
  vars <- seq(p)
  corr <- conf$corr
  B <- Bs[[as.character(p)]]
  nr_targets <- conf$nr_targets
  nr_interventions <- conf$nr_interventions
  method <- as.character(conf$methods)
  effect_size <- conf$effect_sizes

  ints <- list()

  # For each intervention, randomly select the variables to intervene on
  for (j in seq(nr_interventions)) {
    target <- sample(vars, conf$nr_targets)
    int <- rep(0, p)
    int[target] <- 1

    ints[[j]] <- int
  }

  ints[[j + 1]] <- rep(0, p) # observational setting

  # Create error covariance and data
  Sigma <- create_Sigma(p, noise_sd = 1, corr = corr)
  backshift_dat <- create_interventions_data(
    B, Sigma, n, ints, intervention_sd = 1, effect_size = effect_size
  )

  # Estimate causal model using backshift
  warning <- FALSE
  estimate <- tryCatch(
    BACKSHIFT(backshift_dat, n, ev = 0),
    warning = function(warning) {
      warning <<- TRUE
  })

  estimate <- BACKSHIFT(backshift_dat, n, ev = 0)
  estimate[['AhatAdjacency']] <- NULL
  estimate[['percentageOfRunsConverged']] <- NULL

  list(
    'config' = c(time, n, p, corr, nr_interventions, nr_targets, effect_size, warning),
    'estimate' = estimate,
    'Sigmahat' = Sigma
  )
}

end <- Sys.time()
print(end - start)
saveRDS(res, 'Results/backshift-estimates.RDS')
