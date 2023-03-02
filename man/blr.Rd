\name{blr}
\alias{blr}
\title{Bayesian Probit Regression}
\description{
This function implements Bayesian inference for a probit regression model. The function takes as inputs:
- y: a binary response variable
- X: a design matrix of predictors
- offset: an optional offset to the linear predictor
- theta: a vector of initial values for the regression coefficients
- theta_0: a prior mean for the regression coefficients
- N_sim: the number of posterior samples to draw from the model

The function uses a Gibbs sampling algorithm to draw posterior samples from the model. It computes the posterior distribution of the regression coefficients (theta) given the data (y and X) and the prior distribution of the coefficients (theta_0).

The function returns a matrix of N_sim posterior samples for the regression coefficients.
}
\usage{
blr(y, X, offset = 0, theta, theta_0 = 0, N_sim = 1)
}
\arguments{
\item{y}{A binary response variable.}
\item{X}{A design matrix of predictors.}
\item{offset}{An optional offset to the linear predictor.}
\item{theta}{A vector of initial values for the regression coefficients.}
\item{theta_0}{A prior mean for the regression coefficients.}
\item{N_sim}{The number of posterior samples to draw from the model.}
}
\value{
A matrix of N_sim posterior samples for the regression coefficients.
}
\examples{
# generate data
set.seed(1)
n <- 100
X <- matrix(rnorm(n * 2), ncol = 2)
theta_true <- c(1, 2)
z <- X %*% theta_true + rnorm(n)
y <- ifelse(z > 0, 1, 0)

# fit the model
theta_hat <- blr(y, X, theta = c(0,0), N_sim = 1000)

# summary of posterior samples
mean(theta_hat[,1])
mean(theta_hat[,2])
}
\seealso{
\code{\link{rtruncnorm}}
}
\author{
Veronica Vinciotti, Ernst C. Wit and Francisco Richter
}
