bpr <- function(y, X, offset = 0, theta, theta_0 = c(0, 0, 0), N_sim = 1) {
  
  
  # Dimensions of theta
  D <- ncol(X)
  
  # Number of observations
  n <- length(y)
  N1 <- sum(y)
  N0 <- n - N1
  
  # Conjugate prior on the coefficients theta ~ N(theta_0, Q_0)
  Q_0 <- diag(10, D)
  
  # Initialize parameters
  z <- rep(NA, n)
  
  # Matrix storing samples of the theta parameter
  theta_chain <- matrix(0, nrow = N_sim, ncol = D)
  
  # ---------------------------------
  # Gibbs sampling algorithm
  # ---------------------------------
  
  # Compute posterior variance of theta
  prec_0 <- ginv(Q_0)
  V <- ginv(prec_0 + t(X) %*% X)
  
  for (t in 1:N_sim) {
    # Update Mean of z
    mu_z <- X %*% theta + offset
    # Draw latent variable z from its full conditional: z | theta, y, X
    if (N0 > 0) {
      z[y == 0] <- rtruncnorm(N0, mean = mu_z[y == 0], sd = 1, a = -Inf, b = 0)
    }
    if (N1 > 0) {
      z[y == 1] <- rtruncnorm(N1, mean = mu_z[y == 1], sd = 1, a = 0, b = Inf)
    }
    
    # Compute posterior mean of theta
    M <- V %*% (prec_0 %*% theta_0 + t(X) %*% (z - offset))
    # Draw variable theta from its full conditional: theta | z, X
    theta <- rmvnorm(1, M, V)
    
    # Store the theta draws
    theta_chain[t, ] <- theta
  }
  
  return(theta_chain)
}


blr2 <- function(y, X, offset = 0, theta = rep(0, ncol(X)), theta_0 = rep(0, ncol(X)), N_sim = 1) {
  # dim of theta
  D <- ncol(X)
  
  # number of observations
  n <- length(y)
  N1 <- sum(y)
  N0 <- n - N1
  
  # Conjugate prior on the coefficients \theta ~ N(theta_0, Q_0)
  Q_0 <- diag(10, D)
  
  # Initialize parameters
  z <- rep(NA, n)
  
  # Matrix storing samples of the \theta parameter
  theta_chain <- matrix(0, nrow = N_sim, ncol = D)
  
  # ---------------------------------
  # Gibbs sampling algorithm
  # ---------------------------------
  
  # Compute posterior variance of theta
  prec_0 <- solve(Q_0)
  V <- solve(prec_0 + crossprod(X, X))
  
  for (t in 1:N_sim) {
    # Update Mean of z
    mu_z <- X %*% theta + offset
    # Draw latent variable z from its full conditional: z | \theta, y, X
    if (sum(1 - y) > 0) {
      z[y == 0] <- rtruncnorm(N0, mean = mu_z[y == 0], sd = 1, a = -Inf, b = 0)
    }
    if (sum(y) > 0) {
      z[y == 1] <- rtruncnorm(N1, mean = mu_z[y == 1], sd = 1, a = 0, b = Inf)
    }
    
    # Compute posterior mean of theta
    M <- V %*% (prec_0 %*% theta_0 + crossprod(X, z - offset))
    # Draw variable \theta from its full conditional: \theta | z, X
    theta <- c(rmvnorm(1, M, V))
    
    # Store the \theta draws
    theta_chain[t, ] <- theta
  }
  return(theta_chain)
}
