library('dplyr')
library('tidyr')
library('ggplot2')
library('latex2exp')
library('RColorBrewer')
source('Code/helpers.R')
cols <- brewer.pal(3, 'Set1')

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
    d1, d2, d3, d4, ylim, y_axis = TRUE,
    cols = brewer.pal(4, 'Dark2'), ...
  ) {
  
  plot(
    n, d1$mean, axes = FALSE,
    pch = 20, ylim = ylim, cex = 1.25,
    type = 'l', lwd = 2, col = cols[1],
    cex.lab = 1.25, font.main = 1,
    xlab = 'Sample size n', ...
  )
  
  add_sdband(n, d1$mean, d1$sd, cols[1])
  lines(n, d2$mean, lwd = 2, col = cols[2])
  add_sdband(n, d2$mean, d2$sd, col = cols[2])
  
  lines(n, d3$mean, lwd = 2, col = cols[3])
  add_sdband(n, d3$mean, d3$sd, col = cols[3])
  
  lines(n, d4$mean, lwd = 2, col = cols[4])
  add_sdband(n, d4$mean, d4$sd, col = cols[4])
  
  if (y_axis) {
    axis(2, las = 2, at = seq(ylim[1], ylim[length(ylim)], 0.20))
  }
  axis(1, at = c(50, seq(200, 1000, 200)))
}


pdf('Figures/Simulated-Results-Sanity.pdf', width = 9, height = 5)
par(mfrow = c(1, 2))
par(mar = c(5, 5, 4, -0.10) + 0.10)
d1 <- filter(datb, betas == 'B12')
d2 <- filter(datb, betas == 'B23')
d3 <- filter(datb, betas == 'B32')
d4 <- filter(datb, betas == 'B43')
plot_estimation(
  d1, d2, d3, d4, ylim = c(-0.80, 1), ylab = 'Value',
  main = 'Parameter Estimates'
)
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

par(mar = c(5, 1, 4, 3) + 0.10)

d1 <- filter(datc, causal == 'C12')
d2 <- filter(datc, causal == 'C13')
d3 <- filter(datc, causal == 'C43')
d4 <- filter(datc, causal == 'C42')
plot_estimation(
  d1, d2, d3, d4, ylim = c(-0.80, 1), y_axis = FALSE, ylab = '',
  main = 'Causal Effects Estimates',
  cols = brewer.pal(4, 'Set2')
)
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
##### Equilibrium - Snapshot Weights Results
############################################
res <- readRDS('Results/measurement-results.RDS')
datm <- get_dat(res)

title_12 <- expression('Effect of ' ~ X[1] ~ ' on ' ~ X[2])
title_13 <- expression('Effect of ' ~ X[1] ~ ' on ' ~ X[3])
title_42 <- expression('Effect of ' ~ X[4] ~ ' on ' ~ X[2])
title_43 <- expression('Effect of ' ~ X[4] ~ ' on ' ~ X[3])

ylimb <- c(-0.20, 0.80)
ylimb2 <- c(-0.80, 0.20)
par(mfrow = c(2, 2))

pdf('Figures/Simulated-Results-Measurement.pdf', width = 9, height = 7)
par(mfrow = c(2, 2))

par(mar = c(1, 5, 4, -0.10) + 0.10)
plot_box(
  datm, 'X1', 'X2', ylim = ylimb,
  xlab = '', ylab = 'Causal effect',
  main = title_12, x_axis = FALSE
)
axis(2, las = 2, at = seq(ylimb[1], ylimb[length(ylimb)], 0.20))
# label_x <- seq(1.25, by = 1.50, length.out = 6)
# 
# true_effects1 <- comp_press(t(B), rep(0, 4), target = 1, a = 1)
# true_effects4 <- comp_press(t(B), rep(0, 4), target = 4, a = 1)
# lines(x = c(label_x[1], label_x[6]), y = c(true_effects1[2], true_effects1[2]))

legend(
  'topright',
  legend = c(
    expression(rho[12] ~ ' = -0.25'),
    expression(rho[12] ~ ' = 0.00'),
    expression(rho[12] ~ ' = 0.25')
  ),
  fill = brewer.pal(3, 'Dark2'), bty = 'n'
)

par(mar = c(1, 1, 4, 3) + 0.10)
plot_box(
  datm, 'X1', 'X3', ylim = ylimb,
  xlab = '', ylab = '', main = title_13, x_axis = FALSE
)

par(mar = c(4, 5, 1, -0.10) + 0.10)
plot_box(
  datm, 'X4', 'X2', ylim = ylimb2,
  xlab = expression(sigma[epsilon]^2 / (1 + sigma[epsilon]^2)),
  ylab = 'Causal effect',
  main = title_42, x_axis = TRUE
)
axis(2, las = 2, at = seq(ylimb2[1], ylimb2[length(ylimb2)], 0.20))

par(mar = c(4, 1, 1, 3) + 0.10)
plot_box(
  datm, 'X4', 'X3', ylim = ylimb2,
  xlab = expression(sigma[epsilon]^2 / (1 + sigma[epsilon]^2)),
  ylab = '', main = title_43, x_axis = TRUE
)
dev.off()


###################################
##### Equilibrium - Extreme Results
###################################
edgecols <- brewer.pal(3, 'Set1')[-3]
posCol <- edgecols[2]
negCol <- edgecols[1]
pdf('Figures/Simulated-Example.pdf', width = 8, height = 7)
mat <- matrix(c(1, 2, 3, 4, 4, 4), nrow = 2, byrow = TRUE)
nf <- layout(mat, heights = c(0.45, 0.45, 0.10))
par(mar = c(4.5, 4, 4, 0) + 0.10)

B1 <- get_beta(res0, resid_cor = 0, resid_prop = 0)
B2 <- get_beta(res0, resid_cor = 0, resid_prop = 0.50)
B3 <- get_beta(res3b, resid_cor = 0.50, resid_prop = 0.50)
netplot(t(B1), posCol = posCol, negCol = negCol, edge.label.position = c(0.5,0.425,0.425,0.5))
mtext('Equilibrium Data', line=2, cex=1)

par(mar = c(3.2, 0, 4.1, 0), pty = 's')
netplot(t(B2), posCol = posCol, negCol = negCol)
mtext('Snapshot Data', line=2, cex=1)

par(mar = c(3.2, 0, 4.1, 0), pty = 's')
netplot(t(B3), posCol = posCol, negCol = negCol)
mtext('Snapshot Data', line=2, cex=1)

par(mar = c(1, 4, 2, 4) + 0.1, pty = 'm')

dat_res0 <- filter(get_dat(res0), resid_prop == 0)
dat_res0$type <- 'Equilibrium'

dat_res1 <- filter(get_dat(res0), resid_prop == 0.50)
dat_res1$type <- 'Snapshop (uncorrelated)'

dat_res2 <- filter(get_dat(res3b), resid_prop == 0.50, resid_cor == 0.50)
dat_res2$type <- 'Snapshop (correlated)'
dat_all <- rbind(dat_res0, dat_res1, dat_res2)

boxplot(
  value ~ type + target + variable, data = dat_all, axes = FALSE,
  xlim = c(1, 12), ylim = c(-1, 1.5), col = rev(cols),
  ylab = 'Value', xlab = '', pch = 20,
  pars = list(boxwex = 0.50, staplewex = 0.30, outwex = 0.2),
  cex.lab = 1.25, outline = FALSE
)
axis(2, las = 2, cex = 1.5)
legend(
  'bottomleft',
  legend = c('Equilibrium', 'Snapshot (uncorrelated)', 'Snapshot (correlated)'),
  fill = rev(cols), cex = 1, bty = 'n'
)
legend(
  x = 0.575, y = -0.50,
  legend = c('True value'), seg.len = 0.75,
  lwd = 2, col = 'gray', bty = 'n', cex = 1
)
mtext('True and estimated effect of interventions', side = 3, line = 2, cex = 1)
mtext('Setting', side = 1, line = 0, cex = 0.90)

len <- 0.70 / 2
true_values <- c(
  comp_press(t(B), rep(0, 4), target = 1, a = 1)[2:3],
  comp_press(t(B), rep(0, 4), target = 4, a = 1)[2:3]
)
  
lines(c(1 - len, 3 + len), c(true_values[1], true_values[1]), lty = 1, lwd = 2, col = 'gray')
lines(c(4 - len, 6 + len), c(true_values[2], true_values[2]), lty = 1, lwd = 2, col = 'gray')
lines(c(7 - len, 9 + len), c(true_values[3], true_values[3]), lty = 1, lwd = 2, col = 'gray')
lines(c(10 - len, 12 + len), c(true_values[4], true_values[4]), lty = 1, lwd = 2, col = 'gray')

yh <- 1.4
tcex <- 1
text(3.5, yh + 0.10, expression('Targeting ' ~ X[1]), cex = tcex)
text(9.5, yh + 0.10, expression('Targeting ' ~ X[4]), cex = tcex)

yh2 <- yh - 0.30
text(2.0, yh2 + 0.10, expression('Effect on ' ~ X[2]), cex = tcex)
text(5.0, yh2 + 0.10, expression('Effect on ' ~ X[3]), cex = tcex)
text(8.0, yh2 + 0.10, expression('Effect on ' ~ X[2]), cex = tcex)
text(11.0, yh2 + 0.10, expression('Effect on ' ~ X[3]), cex = tcex)
dev.off()


dat$type <- ifelse(dat$resid_prop == 0, 'eq', 'snap')
ratios <- c()
resid_props <- unique(dat$resid_prop)
eqval <- dat[dat$resid_prop == 0, ]$value

for (i in seq(length(unique(resid_props)))) {
  resid_prop <- resid_props[i]
  sel <- dat[dat$resid_prop == resid_prop, ]
  
  ratios <- c(ratios, sel$value / eqval)
}

dat$ratio <- ratios
datsum <- dat %>% 
  filter(resid_prop > 0) %>% 
  group_by(target, variable, resid_prop, resid_cor) %>% 
  summarize(
    mean_ratio = mean(ratio),
    sd_ratio = sd(ratio)
  )

plot_ratio <- function(datsum, target, variable, ylim, main) {
  d <- datsum[datsum$target == target & datsum$variable == variable, ]
  plot(
    d$resid_prop, d$mean_ratio, axes = FALSE,
    pch = 20, ylim = ylim, cex = 1.25,
    cex.lab = 1.25, main = main, xlab = 'Residual proportion', ylab = 'Snapshot / Equilibrium'
  )
  ylo <- d$mean_ratio - d$sd_ratio
  yhi <- d$mean_ratio + d$sd_ratio
  arrows(d$resid_prop, ylo, d$resid_prop, yhi, col = 'black', lwd = 1, length = 0.10, code = 3, angle = 90)
  axis(1, d$resid_prop)
  
  if (!is.null(ylim)) {
    axis(2, las = 2, at = seq(ylim[1], ylim[length(ylim)], 0.10))
  }
}


ylim <- c(-0.20, 1.40)
d <- filter(datsum, resid_cor == 0.00)
par(mfrow = c(2, 2))
plot_ratio(d, 'X1', 'X2', ylim = ylim, title_12)
plot_ratio(d, 'X1', 'X3', ylim = ylim, title_13)
plot_ratio(d, 'X4', 'X2', ylim = ylim, title_42)
plot_ratio(d, 'X4', 'X3', ylim = ylim, title_43)



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
