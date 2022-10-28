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


# Get beta estimates from the simulation results object
get_beta <- function(res, resid_cor = 0, resid_prop = 0) {
  p <- 4
  
  count <- 0
  Best <- matrix(0, p, p)
  
  is_equal <- function(x, y) abs(x - y) < 0.01
  
  for (i in seq(nrow(res))) {
    el <- res[i, ]
    
    if (is_equal(el$resid_cor, resid_cor) && is_equal(el$resid_prop, resid_prop)) {
      count <- count + 1
      Best <- Best + el$Bsnap
    }
  }
  
  Best <- Best / count
  Best
}

# Get results data from the simulation results object
get_dat <- function(res, key = 'dat') {
  dat <- do.call('rbind', lapply(seq(nrow(res)), function(i) {
    res[i, ][[key]]
  }))
  
  dat
}


# Function for creating boxplots
plot_box <- function(dat, target, variable, x_axis = TRUE, ...) {
  cols <- colorRampPalette(brewer.pal(9, 'Blues'))(11)
  cols <- brewer.pal(3, 'Dark2')
  d <- dat[dat$target == target & dat$variable == variable, ]
  
  label_x <- seq(1.25, by = 1.50, length.out = 6)
  at_x <- c(sapply(label_x, function(x) { c(x - 0.40, x, x + 0.40) }))
  
  boxplot(
    value ~ resid_cor + resid_prop, axes = FALSE, data = d,
    col = cols, pch = 20,
    pars = list(boxwex = 0.40, staplewex = 0.40, outwex = 0.20),
    cex.lab = 1.25, cex.main = 1.50, outline = FALSE, at = at_x, ...
  )
  
  if (x_axis) {
    axis(1, at = label_x, labels = unique(dat$resid_prop))
  }
}


# add true effect to plot
add_true_effect <- function(target, effect, xlims = c(50, 1000)) {
  true_effects <- comp_press(B, rep(0, 4), target = target, a = 1)
  lines(
    x = xlims, y = c(true_effects[effect], true_effects[effect]),
    lwd = 2, lty = 2, col = 'gray76'
  )
}

# add true effects for figure 4
add_true_effect2 <- function(target, effect) {
  label_x <- seq(1, by = 1.75, length.out = 6)
  true_effects <- comp_press(B, rep(0, 4), target = target, a = 1)
  lines(
    x = c(0.75, 9.50), y = c(true_effects[effect], true_effects[effect]),
    lwd = 2, lty = 2, col = 'gray76'
  )
}

# add true parameter to plot
add_true_param <- function(B, target, effect, xlims = c(50, 1000)) {
  lines(
    x = xlims, y = c(B[target, effect], B[target, effect]),
    lwd = 2, lty = 2, col = 'gray76'
  )
}

# add standard deviations to plot
add_sdband <- function(n, m, sd, col) {
  polygon(
    c(n, rev(n)), c(m - sd, rev(m + sd)),
    col = adjustcolor(col, 0.25), border = NA
  )
}


# Function for creating line plots
plot_estimation <- function(
    x, d1, d2, d3, d4, ylim, y_axis = TRUE,
    cols = brewer.pal(4, 'Dark2'), ...
) {
  
  plot(
    x, d1$mean, axes = FALSE,
    pch = 20, ylim = ylim, cex = 1.25,
    type = 'l', lwd = 2, col = cols[1],
    cex.lab = 1.25, font.main = 1, ...
  )
  
  add_sdband(x, d1$mean, d1$sd, cols[1])
  lines(x, d2$mean, lwd = 2, col = cols[2])
  add_sdband(x, d2$mean, d2$sd, col = cols[2])
  
  lines(x, d3$mean, lwd = 2, col = cols[3])
  add_sdband(x, d3$mean, d3$sd, col = cols[3])
  
  lines(x, d4$mean, lwd = 2, col = cols[4])
  add_sdband(x, d4$mean, d4$sd, col = cols[4])
  
  if (y_axis) {
    axis(2, las = 2, at = seq(ylim[1], ylim[length(ylim)], 0.20))
  }
}


# version 2 of the compute press function, used in plotting Figure 1
compute_press <- function(c, Phi, target, press_value){
  p <- nrow(Phi)
  Phi_star <- Phi[-target,][,-target]
  c_star <- c[-target]
  Psi <- Phi[,target][-target]
  mu <- rep(NA,p) ;   mu[target] <- press_value
  mu[-target] <- solve(diag(p-1) - Phi_star)%*%(c_star + Psi%*%matrix(press_value,1,1))
  mu
}

# function to create simulated VAR data for Figures 1 - 3
simulate_dynamics <- function(c, B, press = NULL, times = 30, noise = T, noiseSD = 1, start = NULL) {
  p <- nrow(B)
  
  x <- matrix(NA, nrow = times, ncol = p)
  if(is.null(start)){ x[1, ] <- solve(diag(p) - B)%*%c } else {x[1,] <- start }
  
  for(t in seq(2, times)){
    if(noise){epsilon <- rnorm(p,0,noiseSD)} else { epsilon <- rep(0,p)}
    
    x[t, ] <- c + B %*% x[t-1, ] + epsilon
    
    if (!is.null(press)) {
      x[t, ] <- sapply(seq(p), function(i) ifelse(!is.na(press[i]), press[i], x[t, i]))
    }
  }
  
  x
}

# Make a nice looking network plot

netplot <- function(mat, greyscale = FALSE, maximum = .5,asize = 12, edge.labels = TRUE,
                    edge.label.cex = 2, fade = FALSE, shape = "circle",
                    labels = c(expression(X[1]),
                               expression(X[2]),
                               expression(X[3]),
                               expression(X[4])),
                    vsize = 20,
                    esize = 12,
                    posCol = NULL,
                    negCol = NULL, negdash = FALSE, ...){
  
  layout <- rbind(c(0,1), 
                  c(1,1), 
                  c(1,0),
                  c(0,0))
  
  if(isTRUE(labels)){ labels = c("X1", "X2", "X3", "X4") }
  
  m_lty <- matrix(1, 4, 4)
  if(isTRUE(negdash)){m_lty[mat<0] <- 2}
  
  if(is.null(posCol)) posCol <- "firebrick2"
  if(is.null(negCol)) negCol <- "blue"
  
  m_col <- matrix(negCol,4,4)
  m_col[mat > 0 ] <- posCol
  if(greyscale){
    qgraph::qgraph(t(mat), 
                   layout = layout,
                   directed = T, 
                   edge.color = "darkgrey",
                   edge.labels = edge.labels,
                   edge.label.cex = edge.label.cex,
                   edge.label.color = "darkgrey",
                   # curved = FALSE, 
                   lty = t(m_lty), 
                   vsize = vsize, 
                   esize = esize,
                   asize= asize,
                   # color = cols,
                   mar = c(8, 10, 8, 8), maximum=maximum,
                   fade = fade,
                   shape = shape,
                   maximum = maximum,
                   labels = labels, ...)
  } else{
    qgraph::qgraph(t(mat),
                   edge.color = t(m_col),
                   layout = layout,
                   directed = T, 
                   edge.labels = edge.labels,
                   edge.label.cex = edge.label.cex,
                   # curved = FALSE, 
                   lty = t(m_lty), 
                   vsize = 20, 
                   esize = 12,
                   asize= asize,
                   # color = cols,
                   mar = c(8, 10, 8, 8), maximum=maximum,
                   fade = fade,
                   shape = shape,
                   maximum = maximum,
                   labels = labels, ...)
  }
}

# -------------------------------------------

# Time series plot


plotTimeSeries <- function(x, xseq = NULL, cols, ylim = c(-3,3), 
                           xlim = NULL, mu = NULL, 
                           mu_mat  = NULL, mu_seq = NULL, axl = NA, cex.axis = 1,
                           alphamu = 1, alphax = 1, mulwd = 2, cex.lab = 1.5){
  if(is.null(xlim)) xlim = c(0,(nrow(x)))
  if(is.null(xseq)) xseq = 0:(nrow(x)-1)
  plot.new()
  plot.window(xlim = xlim, ylim = ylim)
  title(xlab = "Time (t)", ylab = expression(X[t]), cex.lab = cex.lab, line = axl)
  axis(1, cex.axis = cex.axis); axis(2, cex.axis = cex.axis)
  for(i in 1:ncol(x)){
    lines(xseq,x[,i], col = alpha(cols[i],alphax), lwd = 2)
    if(is.null(mu_mat)){
      abline(h = mu[i], col = alpha(cols[i],alphamu), lty = 2, lwd = mulwd)
    } else{
      lines(y = mu_mat[,i], x = mu_seq, col = alpha(cols[i],alphamu), lty = 2, lwd = mulwd)
    }
    # abline(h = mean(x[,i]), col = cols[i], lty = 3, lwd = 2)
  }
}

plotLegend <- function(x = -1.25, y = .75, lty = c(1,1,1,1,2), cols, 
                       legend = c(expression(X["1,t"]),
                                  expression(X["2,t"]),
                                  expression(X["3,t"]),
                                  expression(X["4,t"]),
                                  expression(mu)),
                       cex = 1.15, lwd = 2,plotnew = TRUE, ...){
  if(isTRUE(plotnew)){
    plot.new()
  }
  par(xpd=TRUE)
  legend(x=x,y = y, lty = lty, col = c(cols, "gray"), 
         legend = legend,
         cex = cex, lwd = lwd, ...)
  par(xpd=FALSE)
}


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
