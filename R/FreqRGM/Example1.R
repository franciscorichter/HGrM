 Load required libraries
library(mvtnorm)
library(stats)
library(Matrix)
library(pracma)

# Set seed for reproducibility
set.seed(42)

# Generate synthetic data
n <- 100  # number of observations
p <- 10   # number of variables

# True graph structure: a simple chain graph
true_G <- matrix(0, p, p)
for (i in 1:(p - 1)) {
  true_G[i, i + 1] <- 1
  true_G[i + 1, i] <- 1
}

# Generate random precision matrix K based on true graph structure
K_true <- optimize_K(S = cov(matrix(rnorm(n * p), n, p)), G = true_G)

# Generate multivariate normal data using the true precision matrix
Sigma_true <- solve(K_true)
mu <- rep(0, p)
Z <- rmvnorm(n, mean = mu, sigma = Sigma_true)

# Generate features for the latent probit model
w <- matrix(rnorm(p * p), p, p)

# Set the penalty parameter lambda
lambda <- 1

# Estimate the graph structure and precision matrix using graph_search
result <- graph_search(Z, w, lambda)

# Extract the estimated graph, precision matrix, and parameters
estimated_G <- result$G
estimated_K <- result$K
estimated_theta <- result$theta

# Compare the true and estimated graph structures
cat("True graph structure:\n")
print(true_G)
cat("Estimated graph structure:\n")
print(estimated_G)

# Compare the true and estimated precision matrices
cat("True precision matrix (K_true):\n")
print(K_true)
cat("Estimated precision matrix (estimated_K):\n")
print(estimated_K)

# Print the estimated parameters
cat("Estimated parameters (theta):\n")
print(estimated_theta)
