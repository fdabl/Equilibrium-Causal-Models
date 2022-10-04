library('MASS')
library('dplyr')
library('lavaan')
library('ggplot2')
library('mvtnorm')
library('backShift')
library('gridExtra')
library('matrixcalc')
library('doParallel')
library('RColorBrewer')

# custom ggplot theme
custom_theme <- theme_minimal() +
  theme(
    legend.position = 'none',
    axis.ticks.y = element_blank(),
    panel.border = element_blank(),
    panel.spacing.x = unit(1, 'lines'),
    panel.spacing.y = unit(1, 'lines'),
    strip.text.x = element_text(size = 14),
    strip.text.y = element_text(size = 14),
    legend.text = element_text(size = 10),
    legend.title = element_text(size = 14),
    axis.text = element_text(size = 10),
    axis.title = element_text(size = 12),
    plot.title = element_text(hjust = 0.5, size = 18),
    plot.subtitle = element_text(hjust = 0.50, size = 14),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_line(colour = "black")
  )


# Rescale VAR(1) matrix to result in equilibrium matrix
rs_beta <- function(B){
  p <- ncol(B)
  ind <- which(diag(B)!=0)
  if(length(ind)== 0){ return(B)} else {
    B_t <- B
    
    # Hytinnen - apply the following repeatedly for all B[i,i] != 0
    for(i in ind){
      # let U be a square 0 matrix with U[i,i] = 1 iff B[i,i] != 0
      U <- matrix(0,p,p)
      U[i,i] <- 1
      B_t <- B_t - (B_t[i,i]/(1-B_t[i,i]))*U%*%(diag(p)-B_t)
    }
    return(B_t)
  }  
}


# Create equilibrium data
create_eqdat <- function(B, intercepts) {
  p <- nrow(B)
  data.frame(intercepts %*% solve(diag(p) - B))
}


# Create noise data
create_noisedat <- function(B, intercepts, Sigma_noise) {
  p <- nrow(B)
  
  Sigma_vec <- solve(diag(p^2) - kronecker(B, B)) %*% vec(Sigma_noise)
  Sigma <- matrix(Sigma_vec, p, p)
  
  dat_noise <- t(apply(intercepts, 1, function(x) {
    rmvnorm(1, mean = rep(0, p), sigma = Sigma)
  }))
  
  dat_noise
}


# Compute the effect of a press intervention
comp_press <- function(Beta, intercepts, target, a = 1){
  p <- nrow(Beta)
  Pk <- diag(p)
  Pk[target,target] <- 0
  avec <- rep(0,p)
  avec[target] <- a
  
  if(any(abs(Re(eigen(Pk%*%Beta)$values))>1)) stop("Intervened System Unstable")
  solve(diag(p)  - Pk%*%Beta)%*%(Pk%*%intercepts + avec)
}


# Compute the effect of a shift intervention
comp_int <- function(Beta, 
                     Sigma = NULL, 
                     x_int = NULL, # (input 1) column vector of intervention values
                     iv = NULL, # (input 2) target variable for intervention
                     effectsize = 1, # (input 2) how many sd's to raise the intercept?
                     cvec = NULL # (input 2) optional vector of intercepts of original system
){
  p <- ncol(Beta)
  if(!is.null(x_int)){
    warning("x_int supplied so ignoring iv, cvec, Sigma and effectsize arguments")
    out <- solve(diag(p)  - Beta)%*%x_int
  } else{
    # if you supplied cvec, i interpret iv and effectsize as CHANGING cvec by a certain amount
    if(!is.null(cvec)) x_int <- cvec else x_int <- rep(0,p)
    # interventions are here defined as changing the intercept by effectsize SDs
    x_int[iv] <- x_int[iv] + sqrt(Sigma[iv,iv])*effectsize
    out <- solve(diag(p) - Beta)%*%x_int
  }
  return(out)
}

# Calculate metrics from estimated backshift matrices
get_met <- function(Bhat, Btrue, Sigmahat, Sigma, max_p = 4) {
  p <- ncol(Btrue)
  
  # Continuous metrics
  L1 <- mean(abs(Bhat - Btrue), na.rm = TRUE)
  L2 <- mean((Bhat - Btrue)^2, na.rm = TRUE)
  
  get_mult <- function(Bhat, Btrue) {
    mult <- abs(Bhat / Btrue)
    mult_median <- median(mult, na.rm = TRUE)
    mult_mean <- mean(mult, na.rm = TRUE)
    mult_min <- min(mult, na.rm = TRUE)
    mult_max <- max(mult, na.rm = TRUE)
    
    c(mult_median, mult_mean, mult_min, mult_max)
  }
  
  mult <- get_mult(Bhat, Btrue)
  
  diag(Bhat) <- 0
  diag(Btrue) <- 0
  
  L1_interventions <- rep(NA, max_p)
  Mult_interventions_mean <- rep(NA, max_p)
  Mult_interventions_median <- rep(NA, max_p)
  Mult_interventions_min <- rep(NA, max_p)
  Mult_interventions_max <- rep(NA, max_p)
  
  for (i in seq(p)) {
    int_est <- comp_int(Bhat, Sigmahat, iv = i, effectsize = 1)
    int_true <- comp_int(Btrue, Sigma, iv = i, effectsize = 1)
    
    L1_interventions[i] <- mean(abs(int_true - int_est))
    mult_int <- abs(int_est / int_true)
    
    Mult_interventions_median[i] <- median(mult_int)
    Mult_interventions_mean[i] <- mean(mult_int)
    Mult_interventions_min[i] <- min(mult_int)
    Mult_interventions_max[i] <- max(mult_int)
  }
  
  c(
    L1, L2, mult, L1_interventions, Mult_interventions_median,
    Mult_interventions_mean, Mult_interventions_min, Mult_interventions_max
  )
}


# Generate data with a particular trait / state variance composition
get_measurement_dat <- function(n, trait_prop = 0.70, resid_cor = NULL) {
  
  # Beta matrix elements
  a <- 0.5
  b <- 1/3
  c <- 0.25
  d <- -0.5
  
  p <- 4
  
  # set latent variances
  p11 <- p44 <- trait_prop
  p22 <- (
    -((a^2)*p11 + (a^2)*b*c*p11 - trait_prop +
        (b^2)*trait_prop +b*c*trait_prop - (b^3)*c*trait_prop)/(1 + b*c)
  )
  
  p33 <- (
    -((d^2)*p44 + b*c*(d^2)*p44 - trait_prop +
        b*c*trait_prop +(c^2)*trait_prop - b*(c^3)*trait_prop)/(1 + b*c)
  )
  
  # Set measurement errors:
  # if standardized variables, total variance is one, so
  theta11 <- theta22 <- theta33 <- theta44 <- 1 - trait_prop
  
  # Make relevant matrices
  I <- diag(p)
  Psi <- diag(c(p11, p22, p33, p44))
  Theta <- diag(c(theta11, theta22, theta33, theta44))
  
  # Add residual correlation between X1 and X2
  if (!is.null(resid_cor)) {
    Psi[1, 2] <- Psi[2, 1] <- resid_cor * sqrt(p11) * sqrt(p22)
  }
  
  # SEM expression for covariance matrix
  S <- (solve(I - B)) %*% Psi %*% t(solve(I - B)) + Theta
  
  dat <- mvrnorm(n, rep(0, p), S)
  colnames(dat) <- paste0('X', seq(p))
  dat
}


# Fit a lavaan model using knowledge about the latent state / trait variance composition
get_lavaan_measurement <- function(dat, trait_prop) {
  
  state_prop <- 1 - trait_prop
  model_spec <- str_glue(
    '
    # latent variable definitions 
    l1 =~ X1
    l2 =~ X2
    l3 =~ X3
    l4 =~ X4
    
    # regressions
    l2 ~ a*l1 + b*l3
    l3 ~ c*l2 + d*l4
    
    # constraints (set ME to zero)
    X1 ~~ {state_prop}*X1
    X2 ~~ {state_prop}*X2
    X3 ~~ {state_prop}*X3
    X4 ~~ {state_prop}*X4
    l1 ~~ {trait_prop}*l1
    l2 ~~ ResVar2*l2
    l3 ~~ ResVar3*l3
    l4 ~~ {trait_prop}*l4
    
    # residual variance specification
    ResVar2 == (
      -((a^2)*{trait_prop} + (a^2)*b*c*{trait_prop} -
      {trait_prop} + (b^2)*{trait_prop} + b*c*{trait_prop} - (b^3)*c*{trait_prop})/(1 + b*c)
    )
    
    ResVar3 == (
    -((d^2)*{trait_prop} + b*c*(d^2)*{trait_prop} -
    {trait_prop} + b*c*{trait_prop} + (c^2)*{trait_prop} - b*(c^3)*{trait_prop})/(1 + b*c)
    )'
  )
  
  fit <- sem(model_spec, data = dat, std.ov = FALSE)
  Bhat <- lavInspect(fit, what = 'est')$beta
  Sigmahat <- lavInspect(fit, what = 'est')$psi
  class(Sigmahat) <- 'matrix'
  
  list('Bhat' = Bhat, 'Sigmahat' = Sigmahat)
}

# Create lavaan model structure from B matrix
create_model_structure <- function(B, fix_coefs = FALSE, estimate_sigma = FALSE) {
  p <- ncol(B)
  ms <- ''
  
  for (i in seq(p)) {
    # Variable i is being influenced by all non-zero elements in this row
    row <- B[i, ]
    
    no_parents <- all(row == 0)
    parents_ix <- which(row != 0)
    
    if (no_parents) {
      string <- paste0('X', i, ' ~ ', '1')
      
    } else {
      
      if (fix_coefs) {
        string <- paste0(
          'X', i, ' ~ 1 + ',
          paste0(B[i, parents_ix], '*X', parents_ix, collapse = ' + ')
        )
        
      } else {
        string <- paste0('X', i, ' ~ 1 + ', paste0('X', parents_ix, collapse = ' + '))
      }
    }
    
    ms <- paste0(ms, string, '\n')
  }
  
  if (estimate_sigma) {
    for (i in seq(p)) {
      for (j in seq(p)) {
        if (i  > j) {
          ms <- paste0(ms, paste0('X', i, ' ~~ ', 'X', j), '\n')
        }
      }
    }
  }
  
  
  ms
}


# Get estimated matrices from lavaan
get_lavaan <- function(dat, B, fix_coefs = FALSE, estimate_sigma = FALSE) {
  cnames <- colnames(dat)
  ms <- create_model_structure(B, fix_coefs = fix_coefs, estimate_sigma = estimate_sigma)
  fit <- sem(ms, std.ov = FALSE, std.lv = FALSE, data = dat)
  Bhat <- lavInspect(fit, what = 'est')$beta[cnames, cnames]
  Sigmahat <- lavInspect(fit, what = 'est')$psi[cnames, cnames]
  inthat <- lavInspect(fit, what = 'est')$alpha[cnames, ]
  class(Sigmahat) <- 'matrix'
  
  list('Bhat' = Bhat, 'Sigmahat' = Sigmahat, 'intercepthat' = inthat)
}


# Estimate causal graph and weights using backshift
BACKSHIFT <- function(
    dat_all, n, threshold = 0.75, eps = 0.05,
    full = FALSE, ev = 1, nsims = 100, ...
  ) {
  p <- ncol(dat_all)
  nr_environments <- nrow(dat_all) / n
  indicator <- rep(seq(nr_environments), each = n)
  
  res <- backShift(dat_all, indicator, threshold = threshold, ev = ev, nsim = nsims, ...)
  res
}


# Creates error covariance matrix
create_Sigma <- function(p, noise_sd = 1, corr = 0) {
  
  Sigma <- diag(noise_sd, p)
  Sigma[upper.tri(Sigma)] <- corr * sample(c(-1, 1), p * (p - 1) / 2, replace = TRUE)
  Sigma[lower.tri(Sigma)] <- t(Sigma)[lower.tri(Sigma)]
  
  while (!is.positive.definite(Sigma)) {
    Sigma <- diag(noise_sd, p)
    Sigma[upper.tri(Sigma)] <- corr * sample(c(-1, 1), p * (p - 1) / 2, replace = TRUE)
    Sigma[lower.tri(Sigma)] <- t(Sigma)[lower.tri(Sigma)]
  }
  
  Sigma
}

# Format data to be in LCC format
format_lcc_dat <- function(dat_list, intervention_list) {
  len <- length(dat_list)
  
  llc <- vector(mode = 'list', len)
  
  for (i in seq(len)) {
    d <- dat_list[[i]]
    llc[[i]]$N <- nrow(d)
    llc[[i]]$data <- d
    llc[[i]]$Cx <- cov(d)
    llc[[i]]$e <- ifelse(is.na(intervention_list[[i]]), 0, 1)
  }
  
  llc
}

#' Creates interventional + observational data sets
create_interventions_data <- function(
  B, Sigma, n, interventions,
  intervention_sd = 1, effect_size = 1
) {
  
  p <- ncol(B)
  I <- diag(p)
  dat <- list()
  len <- length(interventions)
  
  for (i in seq(len)) {
    int <- interventions[[i]]
    targets <- which(int != 0)
    nr_targets <- length(targets)
    # intercept <- matrix(rnorm(n * p, 0, sd = noise_sd), nrow = n)
    intercepts <- rmvnorm(n, mean = rep(0, p), sigma = Sigma)
    
    x <- matrix(rnorm(n * p, 0, intervention_sd), nrow = n)
    
    # For each intervention, the variances need to be different
    # Hence we multiply with a draw from an exponential distribution
    multiplier <- rexp(p, rate = 1 / effect_size) * int
    shift <- sweep(x, 2, multiplier, FUN = '*')
    d <- (intercepts + shift) %*% solve(I - B)
    
    dat[[i]] <- d
  }
  
  do.call('rbind', dat)
}
