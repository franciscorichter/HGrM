\name{graph_search}
\alias{graph_search}
\title{Graph search algorithm for estimating optimal graph structure and precision matrix}
\description{
  The \code{graph_search} function estimates the optimal graph structure and precision matrix for a high-dimensional Gaussian graphical model using a graph search algorithm with logistic regression for the conditional independence structure.
}
\usage{
  graph_search(Z, w, lambda)
}
\arguments{
  \item{Z}{a data matrix where each row corresponds to a sample and each column corresponds to a variable.}

  \item{w}{a vector of weights used in the logistic regression for the conditional independence structure.}

  \item{lambda}{a regularization parameter used to control the number of non-zero edges in the estimated graph structure.}
}
\value{
  A list containing the following components:
  \item{G}{the estimated optimal graph structure.}
  \item{K}{the corresponding precision matrix for the estimated graph structure.}
  \item{theta}{the logistic regression coefficients for the estimated conditional independence structure.}
}
\examples{
  library(Matrix)
  library(MASS)

  # Generate random data
  set.seed(123)
  n <- 100
  p <- 10
  Z <- matrix(rnorm(n * p), n, p)

  # Simulate graph structure
  G <- diag(p)
  G[lower.tri(G)] <- rbinom((p * (p - 1)) / 2, 1, 0.5)
  G <- G + t(G)

  # Ensure that G is positive definite
  G <- nearPD(G)$mat

  # Simulate data from high-dimensional Gaussian graphical model
  S <- solve(G)
  Y <- mvrnorm(n, rep(0, p), S)

  # Estimate optimal graph structure and precision matrix
  w <- rep(1, ncol(Z))
  lambda <- 0.1
  result <- graph_search(Y, w, lambda)
}

\references{
  Meinshausen, N., & Buhlmann, P. (2006). High-dimensional graphs and variable selection with the Lasso. \emph{The Annals of Statistics}, 34(3), 1436-1462.
}
