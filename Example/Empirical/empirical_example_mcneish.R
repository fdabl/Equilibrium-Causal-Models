# code to reproduce the empirical example in Appendix D

library(MplusAutomation)
library(readr)
source("ECM.r")
library(lavaan)

# load estimated cross-lagged mediation model from Mplus
# Mplus model estimated using version 8.1
# data available from: https://osf.io/yk3je/
# in particular code assumes you have the objects "Example Data2.csv" and "Example Data2 (no labels).csv"
# mplus code and output files supplied

obj <- readModels(target = "Example/EmpiricalExample/stationary mediation lag1.out")

posterior_draws <- obj$bparameters$factor.score.imputations
# posterior_draws <- rbind(obj$bparameters$valid_draw$`1`, obj$bparameters$valid_draw$`2`)

b_11_ind <- which(grepl("ADHERE.ON.ADHERE&1", colnames(posterior_draws)))
b_22_ind <- which(grepl("PERS.ON.PERS&1", colnames(posterior_draws)))
b_33_ind <- which(grepl("STEPS.ON.STEPS&1", colnames(posterior_draws)))
b_21_ind <- which(grepl("PERS.ON.ADHERE&1", colnames(posterior_draws)))
b_31_ind <- which(grepl("STEPS.ON.ADHERE&1", colnames(posterior_draws)))
b_32_ind <- which(grepl("STEPS.ON.PERS&1", colnames(posterior_draws)))

# populate array of matrices
bsamps <- array(0, dim = c(3,3,nrow(posterior_draws)))
bsamps[1,1,] <- posterior_draws[,b_11_ind]
bsamps[2,2,] <- posterior_draws[,b_22_ind]
bsamps[3,3,] <- posterior_draws[,b_33_ind]
bsamps[2,1,] <- posterior_draws[,b_21_ind]
bsamps[3,1,] <- posterior_draws[,b_31_ind]
bsamps[3,2,] <- posterior_draws[,b_32_ind]

# obtain point estimates
beta_postmean <- apply(bsamps,c(1,2), mean)

cpars <- obj$parameters$unstandardized
b_11<- cpars[1,"est"]
b_22 <- cpars[2,"est"]
b_21 <- cpars[3,"est"]
b_33 <- cpars[4,"est"]
b_32 <- cpars[5, "est"]
b_31 <- cpars[6,"est"]

# matrix of cross-lagged parameters
beta_pointest <- matrix(c(b_11, 0,0,
                          b_21, b_22, 0,
                          b_31, b_32, b_33), 3,3,byrow = TRUE)


# can also obtain lower and upper bbounds in the same way
b_11_lower <- cpars[1,"lower_2.5ci"]; b_11_upper <- cpars[1,"upper_2.5ci"]
b_22_lower <- cpars[2,"lower_2.5ci"]; b_22_upper <- cpars[2,"upper_2.5ci"]
b_21_lower <- cpars[3,"lower_2.5ci"]; b_21_upper <- cpars[3,"upper_2.5ci"]
b_33_lower <- cpars[4,"lower_2.5ci"]; b_33_upper <- cpars[4,"upper_2.5ci"]
b_32_lower <- cpars[5,"lower_2.5ci"]; b_32_upper <- cpars[5,"upper_2.5ci"]
b_31_lower <- cpars[6,"lower_2.5ci"]; b_31_upper <- cpars[6,"upper_2.5ci"]


lag_beta_list <- rbind(c(b_21, b_21_lower, b_21_upper),
                       c(b_31, b_31_lower, b_31_upper),
                       c(b_32, b_32_lower, b_32_upper),
                       c(b_11, b_11_lower, b_11_upper),
                       c(b_22, b_22_lower, b_22_upper),
                       c(b_33, b_33_lower, b_33_upper))
rownames(lag_beta_list) <- c("X -> M", "X -> Y", "M -> Y", "X -> X", "M -> M", "Y -> Y")
colnames(lag_beta_list) <- c("est","ci.lower","ci.upper")   



# ---------- obtain model-implied ECM  ----------------
# from point estimates
implied_ecm <- rs_beta(beta_pointest)

# now obtain CIs for this by computing the implied ecm from MCMC samples
posterior_ecms <- array(0, dim = c(3,3,nrow(posterior_draws)))
for(i in 1:nrow(posterior_draws)){
  posterior_ecms[,,i] <- rs_beta(bsamps[,,i])
}

implied_ecm_list <- rbind(c(implied_ecm[2,1],quantile(posterior_ecms[2,1,], c(0.025, 0.975))),
                          c(implied_ecm[3,1],quantile(posterior_ecms[3,1,], c(0.025, 0.975))),
                          c(implied_ecm[3,2],quantile(posterior_ecms[3,2,], c(0.025, 0.975))))
implied_ecm_list <- round(implied_ecm_list, 3)

rownames(implied_ecm_list) <- c("X -> M", "X -> Y", "M -> Y")
colnames(implied_ecm_list) <- c("est","ci.lower","ci.upper")

#-----------------------------------------------------------
#------------------ Direct ECM Estimation ------------------
#-----------------------------------------------------------


# ---- extract equilibrium data ----
exdata <- read.csv("Example/EmpiricalExample/Example Data2.csv")
# write to numeric, replace . with NA as missingness indicator
exdata <- as.data.frame(apply(exdata, 2, as.numeric))

# loop through participants, taking time-means of each variable
ids <- unique(exdata$ID)
eq_data_raw <- sapply(ids, function(i){
  tmp <- exdata[exdata$ID == i, ]
  c(mean(tmp$X, na.rm = TRUE), mean(tmp$M, na.rm = TRUE), mean(tmp$Y, na.rm = TRUE))
})

eq_data <- t(eq_data_raw)
colnames(eq_data) <- c("X","M","Y")

# -------- fit equilibrium model ------
syntax_eq <- '
M ~  b21*X
Y ~ b31*X + b32*M
'
fit_eq <- sem(syntax_eq, data = eq_data, meanstructure = TRUE)
summary(fit_eq)

ecm_est<- (lavInspect(fit_eq, what = 'est')$beta)
ecm_est <- ecm_est[c("X","M","Y"), c("X","M","Y")]
ecm_est
implied_ecm

# obtain credible intervals around 
fullpars <- parameterestimates(fit_eq)

estimated_ecm_list <- rbind(fullpars[fullpars$label=="b21",c("est","ci.lower","ci.upper")],
                            fullpars[fullpars$label=="b31",c("est","ci.lower","ci.upper")],
                            fullpars[fullpars$label=="b32",c("est","ci.lower","ci.upper")]
)

rownames(estimated_ecm_list) <- c("X -> M", "X -> Y", "M -> Y")
colnames(estimated_ecm_list) <- c("est","ci.lower","ci.upper")

# show alongside
implied <- matrix(sprintf("%.2f", round(implied_ecm_list,2)),3,3,byrow = F)
estimated_ecm_list <- as.matrix(estimated_ecm_list)
estimated <-  matrix(sprintf("%.2f",round(estimated_ecm_list,2)),3,3,byrow = F)

# latex table output
outmat <- cbind(paste0(implied[,1], " (", implied[,2], ",", implied[,3], ")"),
                paste0(estimated[,1], " (", estimated[,2], ",", estimated[,3], ")")
)

row_names <- c("X -> M", "X -> Y", "M -> Y")
col_names <- c("Implied ECM", "Directly Estimated ECM")

latex_table <- xtable::xtable(outmat,
                              caption = "Point estimates and 95% confidence intervals in parenthesis for the ECM model
                              implied by estimated multilevel VAR(1) model (left column) and directly estimated from
                              the estimated person-specific equilibrium positions.",
                              label = "tab:ECMempirical",
                              row.names = row_names,
                              col.names = col_names)

latex_table

