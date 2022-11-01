library('dplyr')
library('tidyr')
library('qgraph')
library('scales')
library('ggplot2')
library('latex2exp')
library('RColorBrewer')
source('ECM.R')
source('helpers.R')

# set up basic plotting parameters
cols <- brewer.pal(3, 'Set1')
posCol <- cols[2]
negCol <- cols[1]
colsts <- RColorBrewer::brewer.pal(4, 'Dark2')

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


# ------------------------------------------------------------------------------------
# ---------------------------------- Figure 1 ----------------------------------------
# ------------------------------------------------------------------------------------

p <- 4

B <- matrix(c(.5,0,0,0,
              .3,.4,.2,0,
              0,.2,.2,-.4,
              0,0,0,.4),p,p,byrow = T)

intercepts <- c(.5, 1.25, -2, -1.25)
rs_beta(B)

mu <- solve(diag(p) - B) %*% intercepts
mu

set.seed(1234)
x <- simulate_dynamics(intercepts, B, press = NULL, times = 50, noise = T, noiseSD = .3)
xseq <- 0:(nrow(x)-1)


pdf('Figures/Figure1_network_data.pdf', 10*1.2, 4*1.2)
layout(t(c(1,2,3)), widths = c(1,1,.25))
par(xpd = FALSE)
# calls functions in functions_plotting.R
netplot(B, posCol = posCol,negCol = negCol,
        edge.label.position = c(0.5,0.5,0.5,0.425,0.425,0.5,0.5,0.5))
plotTimeSeries(x, xseq, cols = colsts, mu = mu, axl = 2.2, cex.axis = 1.5, cex.lab = 2.2)
plotLegend(cols = colsts, bty = 'n', cex = 2)
dev.off()


# ------------------------------------------------------------------------------------
# ---------------------------------- Figure 2 ----------------------------------------
# ------------------------------------------------------------------------------------

# Figure panels plotted seperately and then combined in Illustrator

# --------------- Pulse Intervention ---------------

par(mfrow = c(1,3))
set.seed(12)
x0 <- simulate_dynamics(intercepts, B, press = NULL, times = 11, noise = T, noiseSD = .3)

target <- 1
pulse <- 0
start <- x0[nrow(x0),]
start[target] <- pulse

set.seed(123)
x1 <- simulate_dynamics(intercepts, B, press= NULL, times = 20, noise = T, noiseSD = .3, start = start)

pulsedata <- rbind(x0,x1)
pulsetime <- c(0:10,10:29)

pdf('Figures/Figure2_pulseint.pdf', 5*1.4, 4*1.4)
plotTimeSeries(pulsedata, xseq = pulsetime, cols = colsts, mu = mu, axl = 2.4, cex.axis = 1.35, cex.lab = 2)
points(x = 10, y = pulse, pch = 23, cex = 2.5, lwd = 2, col = colsts[target], bg = colsts[target])
dev.off()

# --------------- Shift Intervention ---------------

shift_value <- -.5
shift <- intercepts
shift[target] <- intercepts[target] + shift_value
munew <- solve(diag(p) - B)%*%shift
start_shift <- x0[nrow(x0),]


set.seed(123)
xshift <- simulate_dynamics(shift, B, press= NULL, times = 20, noise = TRUE, noiseSD = .3, start = start_shift)
shiftdata <- rbind(x0,xshift)
shifttime <- c(0:10,10:29)
shift_mutime <- c(shifttime,29:30)

mu_mat <- rbind(matrix(rep(mu,11),11,4, byrow = TRUE),
                matrix(rep(munew,22),22,4, byrow = TRUE))

pdf('Figures/Figure2_shift_intervention.pdf',  5*1.4, 4*1.4)
plotTimeSeries(x = shiftdata, xseq = shifttime, cols = colsts, mu_mat = mu_mat, mu_seq = shift_mutime,
               axl = 2.4, cex.axis = 1.35, cex.lab = 2)
points(x = shifttime[11], 
       y = shiftdata[11,target], pch = 21, cex = 2.5, lwd = 2, col = 'black', bg = colsts[target] )
dev.off()

# --------------- Press Intervention ---------------

press_value = -.5

mu_press <- compute_press(c = intercepts, B, target =target, press_value = press_value)
startpress <- x0[nrow(x0),]
startpress[target] <- mu_press[target]

set.seed(12345)
xpress <- simulate_dynamics(intercepts, B, press = c(press_value, NA,NA,NA), times = 20, noise = TRUE, noiseSD =.3,
                            start = startpress)

mu_matpress <- rbind(matrix(rep(mu,11),11,4, byrow = TRUE),
                     matrix(rep(mu_press,22),22,4, byrow = TRUE))

pressdata <- rbind(x0,xpress)
presstime <- c(0:10,10:29)
press_mutime <- c(presstime,29:30)

pdf('Figures/Figure2_pressint.pdf', 5*1.4, 4*1.4)
plotTimeSeries(x = pressdata, xseq = presstime, cols = colsts, mu_mat = mu_matpress, mu_seq = press_mutime,
               axl = 2.4, cex.axis = 1.35, cex.lab = 2)
points(x = presstime[11], 
       y = pressdata[11,target], pch = 23, cex = 2.5, lwd = 2, col = 'black', bg = colsts[target])
dev.off()

# ----------------------- Create Legend ----------------
pdf('Figures/Figure2_legends.pdf',  5*1.4, 4*1.4)
plot.new()
plot.window(c(-1,1), c(-1,1))
legend(x= -1,y = .75, pch = 23, col =  colsts[target], pt.bg = colsts[target], lty = 0,
       legend = c(expression(pulse(X['1,t']))),
       pt.cex = 2, lwd = 2, border = 'white', bty = 'n', cex = 1)
legend(x= -.5,y = .75, pch = 23,  col =  'black', pt.bg = colsts[target], lty = 0,
       legend = c(expression(press(X['1']))),
       pt.cex = 2, lwd = 2, border = 'white', bty = 'n', cex = 1)

legend(x= 0,y = .75, pch = 21, col =  'black', pt.bg = colsts[target], lty = 0,
       legend = c(expression(shift(X['1']))),
       pt.cex = 2, lwd = 2, border = 'white', bty = 'n', cex = 1)
legend(x = .5, y = .75, lty = 2, col = 'gray', legend = expression(mu), cex = 1.15, lwd = 2,
       bty = 'n')

# legend for each line
legend(x = -1.04, y = .5, lty = 1, col = colsts[1], legend = expression(X['1,t']), cex = 1.15, lwd = 2,
       bty = 'n')
legend(x = -.54, y = .5, lty = 1, col = colsts[2], legend = expression(X['2,t']), cex = 1.15, lwd = 2,
       bty = 'n')
legend(x = -0.04, y = .5, lty = 1, col = colsts[3], legend = expression(X['3,t']), cex = 1.15, lwd = 2,
       bty = 'n')
legend(x = .5, y = .5, lty = 1, col = colsts[4], legend = expression(X['4,t']), cex = 1.15, lwd = 2,
       bty = 'n')
dev.off()


# ------------------------------------------------------------------------------------
# ---------------------------------- Figure 3 ----------------------------------------
# ------------------------------------------------------------------------------------
pdf('Figures/Figure3_a.pdf',4,4)
netplot(B, posCol = posCol,negCol = negCol,
        edge.label.position = c(0.5,0.5,0.5,0.425,0.425,0.5,0.5,0.5))
dev.off()

pdf('Figures/Figure3_b.pdf',4,4)
netplot(rs_beta(B), posCol = posCol,negCol = negCol,
        edge.label.position = c(0.5,0.425,0.425,0.5))
dev.off()


# ------------------------------------------------------------------------------------
# ---------------------------------- Figure 5 ----------------------------------------
# ------------------------------------------------------------------------------------

pdf('Figures/Figure5_a.pdf', 5*1.2, 4*1.2)
layout(c(1), widths = c(1))
plotTimeSeries(x, xseq, cols = colsts, mu = mu, axl = 2.4, cex.axis = 1.35, cex.lab = 2,
               alphax = .4, mulwd = 3)
dev.off()

pdf('Figures/Figure5_c.pdf', 5*1.2, 4*1.2)
layout(c(1), widths = c(1))
plotTimeSeries(x, xseq, cols = colsts, mu = mu, axl = 2.4, cex.axis = 1.35, cex.lab = 2,
               alphamu = .5)
t = 45
for(i in 1:4) points(xseq[t], x[t,i], pch = 19, cex = 2, col = alpha(colsts[i], .5), lwd = 2)
dev.off()


# ------------------------------------------------------------------------------------
# -------------------------- Sanity (MLE) Results (Figure 4) -------------------------
# ------------------------------------------------------------------------------------

res <- readRDS('Results/sanity_results.RDS')
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

main <- 'Estimated Parameters Across Sample Size'

pdf('Figures/Figure4_MLE.pdf', width = 9, height = 5)
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

# ------------------------------------------------------------------------------------
# ----------------- Equilibrium - Snapshot Weights Results (Figure 6) ----------------
# ------------------------------------------------------------------------------------

res <- readRDS('Results/measurement_results.RDS')
datm <- get_dat(res)

title_12 <- expression('Effect of ' ~ X[1] ~ ' on ' ~ X[2])
title_13 <- expression('Effect of ' ~ X[1] ~ ' on ' ~ X[3])
title_42 <- expression('Effect of ' ~ X[4] ~ ' on ' ~ X[2])
title_43 <- expression('Effect of ' ~ X[4] ~ ' on ' ~ X[3])

ylimb <- c(-0.20, 0.80)
ylimb2 <- c(-0.80, 0.20)


pdf('Figures/Figure6_Measurement.pdf', width = 9, height = 7)
par(mfrow = c(2, 2))

par(mar = c(1, 5, 4, -0.10) + 0.10)
plot_box(
  datm, 'X1', 'X2', ylim = ylimb,
  xlab = '', ylab = 'Causal effect',
  main = title_12, x_axis = FALSE
)
axis(2, las = 2, at = seq(ylimb[1], ylimb[length(ylimb)], 0.20))
add_true_effect2(1, 2)

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
add_true_effect2(1, 3)

par(mar = c(4, 5, 1, -0.10) + 0.10)
plot_box(
  datm, 'X4', 'X2', ylim = ylimb2,
  xlab = expression(sigma[s]^2 / (1 + sigma[s]^2)),
  ylab = 'Causal effect',
  main = title_42, x_axis = TRUE
)
axis(2, las = 2, at = seq(ylimb2[1], ylimb2[length(ylimb2)], 0.20))
add_true_effect2(4, 2)

par(mar = c(4, 1, 1, 3) + 0.10)
plot_box(
  datm, 'X4', 'X3', ylim = ylimb2,
  xlab = expression(sigma[s]^2 / (1 + sigma[s]^2)),
  ylab = '', main = title_43, x_axis = TRUE
)
add_true_effect2(4, 3)
dev.off()

# ------------------------------------------------------------------------------------
# ------ Equilibrium - Snapshot Weights Latent / State Results (Appendix, Figure 9)  -
# ------------------------------------------------------------------------------------
res <- readRDS('Results/measurement_corr_results.RDS')
datm <- get_dat(res)

pdf('Figures/Figure9_Measurement.pdf', width = 9, height = 7)
par(mfrow = c(2, 2))

par(mar = c(1, 5, 4, -0.10) + 0.10)
plot_box(
  datm, 'X1', 'X2', ylim = ylimb,
  xlab = '', ylab = 'Causal effect',
  main = title_12, x_axis = FALSE
)
axis(2, las = 2, at = seq(ylimb[1], ylimb[length(ylimb)], 0.20))
add_true_effect2(1, 2)

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
add_true_effect2(1, 3)

par(mar = c(4, 5, 1, -0.10) + 0.10)
plot_box(
  datm, 'X4', 'X2', ylim = ylimb2,
  xlab = expression(sigma[s]^2 / (1 + sigma[s]^2)),
  ylab = 'Causal effect',
  main = title_42, x_axis = TRUE
)
axis(2, las = 2, at = seq(ylimb2[1], ylimb2[length(ylimb2)], 0.20))
add_true_effect2(4, 2)

par(mar = c(4, 1, 1, 3) + 0.10)
plot_box(
  datm, 'X4', 'X3', ylim = ylimb2,
  xlab = expression(sigma[s]^2 / (1 + sigma[s]^2)),
  ylab = '', main = title_43, x_axis = TRUE
)
add_true_effect2(4, 3)
dev.off()

# ------------------------------------------------------------------------------------
# ------------------------- State / Trait Results (Figure 7)  ------------------------
# ------------------------------------------------------------------------------------
res <- readRDS('Results/state_trait_results.RDS')
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

pdf('Figures/Figure7_Trait_State.pdf', width = 9, height = 5)
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


# ------------------------------------------------------------------------------------
# --------------- Backshift simulation results (Figures 10 and 11)  ------------------
# ------------------------------------------------------------------------------------
pal <- rev(brewer.pal(9, 'Blues')[-1])
filename <- 'Results/backshift-metrics.csv'

dat <- read.csv(filename) %>% 
  mutate(
    nr_targets_label = factor(nr_targets, labels = paste0('No. targets: ', seq(4))),
    corr_label = factor(corr, labels = c('Confounding: 0', 'Confounding: 0.50')),
    effect_label = factor(effect_size, labels = paste0('Intervention strength: ', unique(effect_size)))
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


# ----------------  Simulation Results Figures ------------

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

pdf('Figures/Figure10_Backshift_Estimation_Error_N.pdf', width = 9, height = 10)
p1
dev.off()

pdf('Figures/Figure11_Backshift_Estimation_Error_Targets.pdf', width = 9, height = 10)
p2
dev.off()
