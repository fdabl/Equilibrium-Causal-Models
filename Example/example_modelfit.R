# Example to show how the measurement model can be 'fixed' in lavaan
# when measurements contain trait and state variance
# --------------------------------------------------------
#  ----------------------- Set-up ------------------------
# --------------------------------------------------------

# load libraries
library('MASS')
library('lavaan')

# load data
x <- readRDS('Example/example_data.RDS')

# ---------------- Lavaan ------------------------

# model syntax, including fixing measurement error variances as described in text
syntax <- '
# latent variable definitions 
l1 =~ x1
l2 =~ x2
l3 =~ x3
l4 =~ x4

# regressions
l2 ~ a*l1 + b*l3
l3 ~ c*l2 + d*l4

# constraints (set ME to zero)
x1 ~~ .3*x1
x2 ~~ .3*x2
x3 ~~ .3*x3
x4 ~~ .3*x4
l1 ~~ .7*l1
l2 ~~ ResVar2*l2
l3 ~~ ResVar3*l3
l4 ~~ .7*l4

# residual variance specification
ResVar2 == -((a^2)*.7+ (a^2)*b*c*.7 - .7 + (b^2)*.7 + b*c*.7 - (b^3)*c*.7)/(1 + b*c)
ResVar3 == -((d^2)*.7 + b*c*(d^2)*.7- .7 + b*c*.7 + (c^2)*.7 - b*(c^3)*.7)/(1 + b*c)
'

# fit the model
fit <- sem(syntax, data = x, std.ov = FALSE)
est <- (lavInspect(fit, what = 'est')$beta)
# approximately the same as the beta matrix used in data generation, see`example_datagen.R`
est
summary(fit)



