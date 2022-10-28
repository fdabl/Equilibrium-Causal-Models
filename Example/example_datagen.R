# file to generate data with "error variances" representing state components 

library(MASS)

# number of variables
q <- 4

# specify beta matrix
a <- 0.5
b<- 0.33
c <- 0.25
d <- -0.5

B <- matrix(c(0,0,0,0,
              a,0,b,0,
              0,c,0,d,
              0,0,0,0),4,4,byrow = T)

# how much variance on the equilibrium/latent level
eqprop <- .7

# set latent variances
# data generated in this way to create a standardized covariance matrix
p11 <- p44 <- eqprop
p22 <- -((a^2)*p11 + (a^2)*b*c*p11 - eqprop + (b^2)*eqprop + b*c*eqprop - (b^3)*c*eqprop)/(1 + b*c)
p33 <- -((d^2)*p44 + b*c*(d^2)*p44 - eqprop + b*c*eqprop + (c^2)*eqprop - b*(c^3)*eqprop)/(1 + b*c)


# set measurement errors variances:
# if standardized variables, total variance is one, so
theta11 <- theta22 <- theta33 <- theta44 <- 1 - eqprop

# make identity matrix
I <- diag(q)

# make other matrices
# Psi; variance covariance matrix of the latent variables
Psi <- diag(c(p11,p22,p33,p44))
# Theta; variance covariance matrix of the measurement errors
Theta <- diag(c(theta11,theta22,theta33,theta44))

# SEM expression for covariance matrix
S = (solve(I - B))%*%Psi%*%t(solve(I-B)) + Theta
S
# diagonal is one, so this worked

# generate data
set.seed(123)
x <- mvrnorm(n = 10000, mu = rep(0,q), Sigma = S)
colnames(x) <- c("x1","x2", "x3","x4")

cov(x)

# save data
saveRDS(x, file = "example_data.RDS")