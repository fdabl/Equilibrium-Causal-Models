# Rescale VAR(1) matrix to equilibrium matrix
rs_beta <- function(B){
  p <- ncol(B)
  ind <- which(diag(B) != 0)
  
  if (length(ind) == 0) {
    return(B)
    
  } else {
    B_t <- B
    
    # Hytinnen - apply the following repeatedly for all B[i,i] != 0
    for (i in ind){
      # let U be a square 0 matrix with U[i,i] = 1 iff B[i,i] != 0
      U <- matrix(0, p, p)
      U[i,i] <- 1
      B_t <- B_t - (B_t[i,i] / (1 - B_t[i,i])) * U %*% (diag(p) - B_t)
    }
    
    return(B_t)
  }  
}


# Compute the effect of a press intervention
comp_press <- function(Beta, intercepts, target, a = 1){
  p <- nrow(Beta)
  Pk <- diag(p)
  Pk[target, target] <- 0
  avec <- rep(0,p)
  avec[target] <- a
  
  if (any(abs(Re(eigen(Pk%*%Beta)$values)) > 1)) {
    stop('Intervened System Unstable')
  }
  
  solve(diag(p)  - Pk %*% Beta) %*% (Pk %*% intercepts + avec)
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
  if (!is.null(x_int)){
    warning('x_int supplied so ignoring iv, cvec, Sigma and effectsize arguments')
    out <- solve(diag(p)  - Beta)%*%x_int
    
  } else{
    # if you supplied cvec, i interpret iv and effectsize as CHANGING cvec by a certain amount
    if (!is.null(cvec)) x_int <- cvec else x_int <- rep(0,p)
    
    # interventions are here defined as changing the intercept by effectsize SDs
    x_int[iv] <- x_int[iv] + sqrt(Sigma[iv,iv]) * effectsize
    out <- solve(diag(p) - Beta) %*% x_int
  }
  
  return(out)
}


# Re-scale the variance covariance matrix of the intercepts
rs_sigma <- function(S, B){
  p <- ncol(S)
  ind <- which(diag(B) ! =0)
  
  if (length(ind) == 0) {
    return(S)
  } else {
    S_t <- S
    # Hytinnen - apply the following repeatedly for all B[i,i] != 0
    for (i in ind){
      # let U be a square 0 matrix with U[i,i] = 1 iff B[i,i] != 0
      U <- matrix(0,p,p) ; I <- diag(p)
      U[i,i] <- 1
      S_t <- (I + (B[i,i] / (1 - B[i,i])) * U) %*% S_t %*% t(I + (B[i,i] / (1 - B[i,i])) * U)
    }
    
    return(S_t)
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
