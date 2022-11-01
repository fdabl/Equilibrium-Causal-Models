library('dplyr')
library('tidyr')
library('ggplot2')
library('latex2exp')
library('RColorBrewer')
source('Code/helpers.R')
cols <- brewer.pal(3, 'Set1')

B4 <- t(matrix(
  c(.5,   0, 0,  0,
    .3, .4, .2,  0,
    0,  .2, .2, -.4,
    0,   0,  0,  .4), 4, 4, byrow = T
))

B <- rs_beta(t(B4))

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



###########################
##### Sanity Results Figure
###########################
res <- readRDS('Results/sanity-results.RDS')
dat <- get_dat(res)

add_true_effect <- function(target, effect, xlims = c(50, 1000)) {
  true_effects <- comp_press(B, rep(0, 4), target = target, a = 1)
  lines(
    x = xlims, y = c(true_effects[effect], true_effects[effect]),
    lwd = 2, lty = 2, col = 'gray76'
  )
}

add_true_param <- function(B, target, effect, xlims = c(50, 1000)) {
  lines(
    x = xlims, y = c(B[target, effect], B[target, effect]),
    lwd = 2, lty = 2, col = 'gray76'
  )
}

datb <- dat %>% 
  group_by(betas, n) %>% 
  summarize(
    mean = mean(beta_values),
    sd = sd(beta_values)
  )

datc <- dat %>% 
  group_by(causal, n) %>% 
  summarize(
    mean = mean(causal_values),
    sd = sd(causal_values)
  )

n <- unique(datb$n)

add_sdband <- function(n, m, sd, col) {
  polygon(
    c(n, rev(n)), c(m - sd, rev(m + sd)),
    col = adjustcolor(col, 0.25), border = NA
  )
}

main <- 'Estimated Parameters Across Sample Size'

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


pdf('Figures/Simulated-Results-Sanity.pdf', width = 9, height = 5)
par(mfrow = c(1, 2))
par(mar = c(5, 5, 4, -0.10) + 0.10)
d1 <- filter(datb, betas == 'B12')
d2 <- filter(datb, betas == 'B23')
d3 <- filter(datb, betas == 'B32')
d4 <- filter(datb, betas == 'B43')
plot_estimation(
  d1$n, d1, d2, d3, d4, ylim = c(-0.80, 1), ylab = 'Value',
  main = 'Parameter Estimates', xlab = 'Sample size n'
)
axis(1, at = c(50, seq(200, 1000, 200)))
add_true_param(t(B), 1, 2)
add_true_param(t(B), 2, 3)
add_true_param(t(B), 3, 2)
add_true_param(t(B), 4, 3)
legend(
  'topleft',
  legend = c(
    expression(beta[12]),
    expression(beta[23]),
    expression(beta[32]),
    expression(beta[43])
  ),
  ncol = 2, x.intersp = 0.50,
  fill = brewer.pal(4, 'Dark2'), bty = 'n'
)
legend(
  'topright',
  legend = c('True value'),
  lwd = 2, lty = 2, col = 'gray76', bty = 'n'
)

par(mar = c(5, 1, 4, 3) + 0.10)

d1 <- filter(datc, causal == 'C12')
d2 <- filter(datc, causal == 'C13')
d3 <- filter(datc, causal == 'C43')
d4 <- filter(datc, causal == 'C42')
plot_estimation(
  d1$n, d1, d2, d3, d4, ylim = c(-0.80, 1), y_axis = FALSE, ylab = '',
  main = 'Causal Effects Estimates',
  xlab = 'Sample size n', cols = brewer.pal(4, 'Set2')
)
axis(1, at = c(50, seq(200, 1000, 200)))
add_true_effect(1, 2)
add_true_effect(1, 3)
add_true_effect(4, 3)
add_true_effect(4, 2)
legend(
  'topleft',
  legend = c(
    expression(X[1] %->% X[2]),
    expression(X[1] %->% X[3]),
    expression(X[4] %->% X[3]),
    expression(X[4] %->% X[2])
  ),
  ncol = 2, x.intersp = 0.50,
  fill = brewer.pal(4, 'Set2'), bty = 'n'
)
dev.off()


############################################
##### Equilibrium - Snapshot Weights Results (Figure 6)
############################################
res <- readRDS('Results/measurement-results.RDS')
datm <- get_dat(res)

title_12 <- expression('Effect of ' ~ X[1] ~ ' on ' ~ X[2])
title_13 <- expression('Effect of ' ~ X[1] ~ ' on ' ~ X[3])
title_42 <- expression('Effect of ' ~ X[4] ~ ' on ' ~ X[2])
title_43 <- expression('Effect of ' ~ X[4] ~ ' on ' ~ X[3])

ylimb <- c(-0.20, 0.80)
ylimb2 <- c(-0.80, 0.20)

add_true_effect <- function(target, effect) {
  label_x <- seq(1, by = 1.75, length.out = 6)
  true_effects <- comp_press(B, rep(0, 4), target = target, a = 1)
  lines(
    x = c(0.75, 9.50), y = c(true_effects[effect], true_effects[effect]),
    lwd = 2, lty = 2, col = 'gray76'
  )
}

pdf('Figures/Simulated-Results-Measurement.pdf', width = 9, height = 7)
par(mfrow = c(2, 2))

par(mar = c(1, 5, 4, -0.10) + 0.10)
plot_box(
  datm, 'X1', 'X2', ylim = ylimb,
  xlab = '', ylab = 'Causal effect',
  main = title_12, x_axis = FALSE
)
axis(2, las = 2, at = seq(ylimb[1], ylimb[length(ylimb)], 0.20))
add_true_effect(1, 2)

legend(
  'bottomleft',
  legend = c(
    expression(rho[12] ~ ' = -0.25'),
    expression(rho[12] ~ ' = 0.00'),
    expression(rho[12] ~ ' = 0.25')
  ),
  fill = brewer.pal(3, 'Dark2'), bty = 'n'
)
legend(
  x = -0.15, y = 0.125,
  legend = c('True effect'),
  lwd = 2, lty = 2, col = 'gray76', bty = 'n'
)

par(mar = c(1, 1, 4, 3) + 0.10)
plot_box(
  datm, 'X1', 'X3', ylim = ylimb,
  xlab = '', ylab = '', main = title_13, x_axis = FALSE
)
add_true_effect(1, 3)

par(mar = c(4, 5, 1, -0.10) + 0.10)
plot_box(
  datm, 'X4', 'X2', ylim = ylimb2,
  xlab = expression(sigma[s]^2 / (1 + sigma[s]^2)),
  ylab = 'Causal effect',
  main = title_42, x_axis = TRUE
)
axis(2, las = 2, at = seq(ylimb2[1], ylimb2[length(ylimb2)], 0.20))
add_true_effect(4, 2)

par(mar = c(4, 1, 1, 3) + 0.10)
plot_box(
  datm, 'X4', 'X3', ylim = ylimb2,
  xlab = expression(sigma[s]^2 / (1 + sigma[s]^2)),
  ylab = '', main = title_43, x_axis = TRUE
)
add_true_effect(4, 3)
dev.off()


######################################################################
##### Equilibrium - Snapshot Weights Latent / State Results (Appendix)
######################################################################
res <- readRDS('Results/measurement-corr-results.RDS')
datm <- get_dat(res)

pdf('Figures/Simulated-Results-Measurement-Latent-State.pdf', width = 9, height = 7)
par(mfrow = c(2, 2))

par(mar = c(1, 5, 4, -0.10) + 0.10)
plot_box(
  datm, 'X1', 'X2', ylim = ylimb,
  xlab = '', ylab = 'Causal effect',
  main = title_12, x_axis = FALSE
)
axis(2, las = 2, at = seq(ylimb[1], ylimb[length(ylimb)], 0.20))
add_true_effect(1, 2)

legend(
  'bottomleft',
  legend = c(
    expression(rho[12] ~ ' = -0.25'),
    expression(rho[12] ~ ' = 0.00'),
    expression(rho[12] ~ ' = 0.25')
  ),
  fill = brewer.pal(3, 'Dark2'), bty = 'n'
)
legend(
  x = -0.15, y = 0.125,
  legend = c('True effect'),
  lwd = 2, lty = 2, col = 'gray76', bty = 'n'
)

par(mar = c(1, 1, 4, 3) + 0.10)
plot_box(
  datm, 'X1', 'X3', ylim = ylimb,
  xlab = '', ylab = '', main = title_13, x_axis = FALSE
)
add_true_effect(1, 3)

par(mar = c(4, 5, 1, -0.10) + 0.10)
plot_box(
  datm, 'X4', 'X2', ylim = ylimb2,
  xlab = expression(sigma[s]^2 / (1 + sigma[s]^2)),
  ylab = 'Causal effect',
  main = title_42, x_axis = TRUE
)
axis(2, las = 2, at = seq(ylimb2[1], ylimb2[length(ylimb2)], 0.20))
add_true_effect(4, 2)

par(mar = c(4, 1, 1, 3) + 0.10)
plot_box(
  datm, 'X4', 'X3', ylim = ylimb2,
  xlab = expression(sigma[s]^2 / (1 + sigma[s]^2)),
  ylab = '', main = title_43, x_axis = TRUE
)
add_true_effect(4, 3)
dev.off()

##################################
##### State / Trait Results Figure
##################################
res <- readRDS('Results/state-trait-results.RDS')
dat <- get_dat(res)

datb <- dat %>% 
  group_by(betas, trait_prop) %>% 
  summarize(
    mean = mean(beta_values),
    sd = sd(beta_values)
  )

datc <- dat %>% 
  group_by(causal, trait_prop) %>% 
  summarize(
    mean = mean(causal_values),
    sd = sd(causal_values)
  )

trait_props <- unique(datb$trait_prop)
main <- 'Estimated Parameters Across Sample Size'

pdf('Figures/Simulated-Results-Trait-State.pdf', width = 9, height = 5)
par(mfrow = c(1, 2))
par(mar = c(5, 5, 4, -0.10) + 0.10)
d1 <- filter(datb, betas == 'B12')
d2 <- filter(datb, betas == 'B23')
d3 <- filter(datb, betas == 'B32')
d4 <- filter(datb, betas == 'B43')
plot_estimation(
  d1$trait_prop, d1, d2, d3, d4, ylim = c(-0.80, 1), ylab = 'Value',
  main = 'Parameter Estimates', xlab = 'Modeled trait variance'
)
xlim <- c(0.50, 0.90)
axis(1, at = seq(0.50, 0.90, 0.10))
#lines(x = c(0.70, 0.70), y = c(-0.80, 1), lwd = 2, col = 'gray76', lty = 2)
lines(x = c(0.70, 0.70), y = c(-0.80, 1), lwd = 1, col = 'black', lty = 1)
add_true_param(t(B), 1, 2, xlim = xlim)
add_true_param(t(B), 2, 3, xlim = xlim)
add_true_param(t(B), 3, 2, xlim = xlim)
add_true_param(t(B), 4, 3, xlim = xlim)
legend(
  'topleft',
  legend = c(
    expression(beta[12]),
    expression(beta[23]),
    expression(beta[32]),
    expression(beta[43])
  ),
  ncol = 2, x.intersp = 0.50,
  fill = brewer.pal(4, 'Dark2'), bty = 'n'
)

legend(
  'topright',
  legend = c('True value'),
  lwd = 2, lty = 2, col = 'gray76', bty = 'n'
)

par(mar = c(5, 1, 4, 3) + 0.10)

d1 <- filter(datc, causal == 'C12')
d2 <- filter(datc, causal == 'C13')
d3 <- filter(datc, causal == 'C43')
d4 <- filter(datc, causal == 'C42')
plot_estimation(
  d1$trait_prop, d1, d2, d3, d4, ylim = c(-0.80, 1), y_axis = FALSE, ylab = '',
  main = 'Causal Effects Estimates',
  xlab = 'Modeled trait variance',
  cols = brewer.pal(4, 'Set2')
)
axis(1, at = seq(0.50, 0.90, 0.10))
#lines(x = c(0.70, 0.70), y = c(-0.80, 1), lwd = 2, col = 'gray76', lty = 2)
add_true_effect(1, 2, xlim = xlim)
add_true_effect(1, 3, xlim = xlim)
add_true_effect(4, 3, xlim = xlim)
add_true_effect(4, 2, xlim = xlim)
lines(x = c(0.70, 0.70), y = c(-0.80, 1), lwd = 1, col = 'black', lty = 1)

legend(
  'topright',
  legend = c(
    expression(X[1] %->% X[2]),
    expression(X[1] %->% X[3]),
    expression(X[4] %->% X[3]),
    expression(X[4] %->% X[2])
  ),
  ncol = 2, x.intersp = 0.50,
  fill = brewer.pal(4, 'Set2'), bty = 'n'
)
dev.off()



##################################
##### Backshift simulation results
##################################
pal <- rev(brewer.pal(9, 'Blues')[-1])
filename <- 'Results/backshift-metrics.csv'

dat <- read.csv(filename) %>% 
  mutate(
    nr_targets_label = factor(dat$nr_targets, labels = paste0('No. targets: ', seq(4))),
    corr_label = factor(dat$corr, labels = c('Confounding: 0', 'Confounding: 0.50')),
    effect_label = factor(dat$effect_size, labels = paste0('Intervention strength: ', unique(dat$effect_size)))
  )

datplot_a <- filter(dat, effect_size == 1, nr_targets == 3, corr %in% c(0, 0.50), n > 100)
datplot_b <- filter(dat, effect_size == 1, corr %in% c(0, 0.50), n == 500)

# Decrease in estimation error as we go from n = 250 to n = 500
dat %>% 
  filter(
    nr_int == 3, n %in% c(250, 500), p == 4, effect_size == 1
  ) %>% 
  group_by(n, corr) %>% 
  summarize(
    L1_mean = mean(L1)
  )

# Decrease in estimation error as we go from s = 3 to s = 4
dat %>% 
  filter(
    nr_int %in% c(3, 4), n == 250, p == 4, effect_size == 1
  ) %>% 
  group_by(nr_int, corr) %>% 
  summarize(
    L1_mean = mean(L1)
  )

# Decrease in estimation error as we go from t = 1 to t = 4 in s = 4
dat %>% 
  filter(
    nr_int %in% c(1, 4), effect_size == 1, n == 500
  ) %>% 
  group_by(corr, nr_targets) %>% 
  summarize(
    L1_mean = mean(L1)
  )

################################
##### Simulation Results Figures
################################
datplot1 <- filter(dat, corr %in% c(0, 0.50))
datplot2 <- filter(dat, corr %in% c(0, 0.50), n == 500)

p1 <- ggplot(datplot1, aes(x = factor(n), y = L1, fill = factor(nr_int))) +
  geom_boxplot(outlier.size = 0.1) +
  facet_wrap(~ effect_label + corr_label, nrow = 3) +
  scale_y_continuous(breaks = seq(0, 0.60, 0.10), limits = c(0, 0.60), expand = c(0, 0)) +
  xlab('Sample size per setting') +
  ylab('Average absolute difference') +
  guides(fill = guide_legend(title = 'Number of settings', nrow = 2)) +
  scale_fill_manual(values = pal) +
  ggtitle(TeX('Overall estimation error across sample size')) +
  custom_theme +
  theme(legend.position = 'top') +
  guides(
    fill = guide_legend(
      title.position = 'top', nrow = 1,
      title = 'Number of settings', title.hjust = 0.50
    )
  )

p2 <- ggplot(datplot2, aes(x = factor(nr_targets), y = L1, fill = factor(nr_int))) +
  geom_boxplot(outlier.size = 0.1) +
  facet_wrap(~ effect_label + corr_label, nrow = 3) +
  scale_y_continuous(breaks = seq(0, 0.60, 0.10), limits = c(0, 0.60), expand = c(0, 0)) +
  xlab('Number of targets') +
  ylab('Average absolute difference') +
  guides(fill = guide_legend(title = 'Number of settings', nrow = 2)) +
  scale_fill_manual(values = pal) +
  ggtitle(TeX('Overall estimation error across number of targets')) +
  custom_theme +
  theme(legend.position = 'top') +
  guides(
    fill = guide_legend(
      title.position = 'top', nrow = 1,
      title = 'Number of settings', title.hjust = 0.50
    )
  )

pdf('Figures/Backshift-Estimation-Error-N.pdf', width = 9, height = 10)
p1
dev.off()

pdf('Figures/Backshift-Estimation-Error-Targets.pdf', width = 9, height = 10)
p2
dev.off()




