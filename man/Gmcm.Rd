\name{Gmcmc}
\alias{Gmcmc}
\title{Random Graph Model}
\description{
  The function implements Bayesian inference of a random graph model for a given input matrix `G` containing graphs across B different environments. 
  The function returns a list of the posterior distributions of environment-specific intercepts `alpha`, latent environment locations `cloc`, and regression coefficients `beta` (if the input matrix `Z` of edge covariates is provided). 
}
\usage{
Gmcmc(G, Z = NULL, n.iter = 1000, alpha = NULL, beta = NULL, cloc = NULL, n.burnin = 500)
}
\arguments{
  \item{G}{A binary matrix representing the graph structures, with dimension n.edges x B.}
  \item{Z}{An optional matrix of covariates with dimension n.edges x ncol(Z).}
  \item{n.iter}{The number of MCMC iterations.}
  \item{alpha}{An optional vector of initial values for the condition-specific intercepts.}
  \item{beta}{An optional vector of initial values for the regression coefficients associated to the covariates in `Z`.}
  \item{cloc}{An optional matrix of initial values for the latent condition locations.}
  \item{n.burnin}{The number of burn-in iterations.}
}
\value{
  A list with components:
  \item{alpha}{The MCMC samples of the condition-specific intercepts.}
  \item{cloc}{The MCMC samples of the latent condition locations.}
  \item{beta}{The MCMC samples of the regression coefficients associated to the covariates in `Z` (if provided).}
}
\examples{
G <- matrix(rbinom(100, 1, 0.5), ncol = 10)
result <- Gmcmc(G, n.iter = 1000, n.burnin = 500)
}
\references{
  V. J. Anderson and J. Robinson (2008). "Bayesian analysis of binary and polychotomous response data". \emph{Journal of the Royal Statistical Society: Series B}, 70(3), pp. 579-599.
}
\author{
  Veronica Vinciotti, Ernst C. Wit and Francisco Richter
}
\seealso{
  \code{\link{blr}}
}
