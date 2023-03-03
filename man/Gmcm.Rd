\name{Gmcmc}
\alias{Gmcmc}
\title{Random Graph Model}
\description{
  The function implements Bayesian inference of a random graph model for a given input matrix `G` containing graphs across B different environments. 
  The function returns a list of the posterior distributions of environment-specific intercepts `alpha`, latent environment locations `loc`, and regression coefficients `theta` (if the input matrix `X` of edge covariates is provided). 
}
\usage{
Gmcmc(G, X = NULL, iter = 1000, alpha = NULL, theta = NULL, loc = NULL, burnin = 0)
}
\arguments{
  \item{G}{A binary matrix representing the graph structures, with dimension number of edges (n.edge) x B.}
  \item{X}{An optional matrix of covariates with dimension n.edge x ncol(X).}
  \item{iter}{The number of MCMC iterations.}
  \item{alpha}{An optional vector of initial values for the condition-specific intercepts. If \code{NULL}, these are drawn from a N(0,1).}
  \item{theta}{An optional vector of initial values for the regression coefficients associated to the covariates in X. If \code{NULL}, these are set to zero.}
  \item{loc}{An optional matrix of initial values for the latent condition locations. If \code{NULL}, these are drawn from a N(0,1).}
  \item{burnin}{The number of burn-in iterations. Default is 0.}
}
\value{
  A list with components:
  \item{alpha}{The MCMC samples of the condition-specific intercepts.}
  \item{loc}{The MCMC samples of the latent condition locations.}
  \item{theta}{The MCMC samples of the regression coefficients associated to the covariates in X (if provided).}
}
\examples{
G <- matrix(rbinom(100, 1, 0.5), ncol = 10)
result <- Gmcmc(G, iter = 1000, burnin = 500)
}
\references{
  V. J. Anderson and J. Robinson (2008). "Bayesian analysis of binary and polychotomous response data". \emph{Journal of the Royal Statistical Society: Series B}, 70(3), pp. 579-599.
}
\author{
  Veronica Vinciotti, Ernst C. Wit and Francisco Richter
}
\seealso{
  \code{\link{bpr}}
}
