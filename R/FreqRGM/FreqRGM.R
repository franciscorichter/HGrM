# Required library
library(mvtnorm)
library(stats)
library(Matrix)
library(pracma)

library(ggm)

# Likelihood for the precision matrix K given the graph G
likelihood_K <- function(K, n, S) {
  log_det_K <- determinant(K, logarithm = TRUE)$modulus
  trace_SK <- sum(S * K)
  C <- n / 2
  log_likelihood_K <- C * (log_det_K - trace_SK)
  return(log_likelihood_K)
}


add_edge <- function(G, i, j) {
  # Add an edge between nodes i and j in graph G
  G[i, j] <- 1
  G[j, i] <- 1
  return(G)
}

remove_edge <- function(G, i, j) {
  # Remove the edge between nodes i and j in graph G
  G[i, j] <- 0
  G[j, i] <- 0
  return(G)
}

delta_loglik_K_edge_remove<-function(K,i,j,n){
  cond.cor<- -K[i,j]/sqrt(K[i,i]*K[j,j])
  delta <- -n*log(1-cond.cor^2)/2
  return(delta)
}


# Likelihood for the graph G
latent_probit_prob <- function(theta,w) {
  # Calculate probability
  Phi_arg <- alpha + crossprod(w, beta)
  probs <- pnorm(Phi_arg)

  return(probs)
}

likelihood_G <- function(G,theta,w) {
  # Extract the parameters from theta
  alpha <- theta$alpha
  beta <- theta$beta
  #c <- theta$c

  # Calculate the log-likelihood for G using the latent probit model
  prob_matrix <- latent_probit_prob(theta, w)
  log_likelihood <- sum(G * log(prob_matrix) + (1 - G) * log(1 - prob_matrix))

  return(log_likelihood)
}




# Log-likelihood of the observed data Z given the precision matrix K
likelihood_Z <- function(K, Z, mu) {
  n <- dim(Z)[1]
  p <- dim(Z)[2]
  log_det_Sigma <- -n * log(det(K)) / 2
  quadratic_form <- 0
  for (i in 1:n) {
    z_i <- Z[i, ]
    diff <- z_i - mu
    quadratic_form <- quadratic_form + t(diff) %*% K %*% diff
  }
  log_likelihood_Z <- -n * p * log(2 * pi) / 2 + log_det_Sigma - quadratic_form / 2
  return(log_likelihood_Z)
}


log_likelihood <- function(K, G, Z, theta, lambda =0) {
  # Extract the parameters from theta
  alpha <- theta$alpha
  beta <- theta$beta
  #c <- theta$c
  w <- theta$w
  nu <- 1
  V <- diag(ncol(Z))
  mu <- apply(Z,2, mean)

  # Calculate the log-likelihood for G using the latent probit model
  log_likelihood_G_value <- likelihood_G(alpha, beta, c, w, G, Z)

  # Calculate the log-likelihood for the precision matrix K given the graph G
  log_likelihood_K_value <- likelihood_K(K, nu, V, S)

  # Calculate the log-likelihood of the observed data Z given the precision matrix K
  # log_likelihood_Z_value <- likelihood_Z(K, Z, mu)

  # L0 penalty
  penalty <- lambda*sum(G!=0)

  # Calculate the total log-likelihood
  total_log_likelihood <- log_likelihood_G_value + log_likelihood_K_value  + penalty

  return(total_log_likelihood)
}



optimize_K <- function(S, G) {
  n <- dim(S)[1]
  I <- diag(n)
  Gamma <- I * (1 - G)
  K <- solve(S - Gamma)
  return(K)
}



