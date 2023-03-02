\name{sim.rgm}
\alias{sim.rgm}
\title{
Simulate Data from a Random Graphical Model with Binary Edges and a Continuous Covariate
}
\usage{
sim.rgm(n = 1000, D = 2, p = 81, B = 10, seed = 123, mcmc_iter = 50, alpha = NULL, theta = NULL, loc = NULL, X = NULL) 
}
\arguments{
\item{n}{The number of observations for each environment. Default is 1000.}
\item{D}{The dimension of the latent space. Default is 2.}
\item{p}{The number of nodes in each graph. Default is 81.}
\item{B}{The number of conditions. Default is 10.}
\item{seed}{The random seed. Default is 123.}
\item{mcmc_iter}{The number of MCMC sampling for the generation of the graphs from the joint random graph distribution. Default is 50.}
\item{alpha}{The true values of the condition-specific intercepts. If \code{NULL}, these are drawn from a N(-2,1) distribution.}
\item{theta}{The true values of the regression coefficients associated to the covariates in X. If \code{NULL}, this is set to 2.5.}
\item{loc}{The true coordinates of the B locations in the latent space. If \code{NULL}, these are drawn from a N(0,0.3) distribution.}
\item{X}{The edge specific covariates. If \code{NULL}, the data for one covariates is drawn from a Uniform(-0.5,0.5) distribution.}
}
\value{
A list with the following elements:
\item{data}{A list of B elements, where each element contains an n-by-p matrix of simulated Gaussian data.}
\item{X}{A dataset of edge covariates, with dimension number of edges per number of covariates.}
\item{loc}{A matrix of the true condition-specific coordinates. Dimension is nrow = B, ncol = D.}
\item{alpha}{A vector of the true condition-specific intercepts. Length is B.}
\item{theta}{A vector of the true regression coefficients associated to the covariates in X.}
\item{G}{A matrix of the true graphs, of dimension number of edges per number of conditions B.}
\item{diagnostic}{The sparsity of the graphs generated across the \code{mcmc_iter} iterations, as a diagnostic tool for convergence.}
}
\description{
This function simulates data from a random graphical model with binary edges and a continuous covariate. The random graph model includes condition-specific intercepts, a 2D latent space model and an edge specific covariate.
}
\examples{
sim_data <- sim.rgm()
}
