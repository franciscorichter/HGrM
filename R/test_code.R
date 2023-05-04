# Load required packages
#install.packages("huge")
library(huge)
#install.packages("MASS")
library(MASS)

# Generate synthetic dataset with known structure
set.seed(42)
n <- 100  # number of samples
p <- 10   # number of variables
Sigma <- matrix(0, p, p)
for (i in 1:p) {
  for (j in 1:p) {
    Sigma[i, j] <- 0.5^(abs(i - j))
  }
}
true_G <- Sigma != 0
rownames(Sigma) <- colnames(Sigma) <- paste("V", 1:p, sep = "")
Z <- mvrnorm(n, mu = rep(0, p), Sigma = Sigma)

# Create a weights vector w
set.seed(1)
w <- rnorm(p * (p - 1) / 2)

# Set lambda value for regularization
lambda <- 0.1

# Run the graph_search function on the synthetic dataset
result <- graph_search(Z, w, lambda)

# Print the resulting graph and precision matrix
print("Estimated Graph (G):")
print(result$G)
print("Estimated Precision Matrix (K):")
print(result$K)
print("Estimated Regression Coefficients (theta):")
print(result$theta)

