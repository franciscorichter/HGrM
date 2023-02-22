\name{sim.rgm}
\alias{sim.rgm}
\title{
Simulate Data from a Latent Space Model with Binary Edges and a Continuous Covariate
}
\usage{
sim.rgm(D = 2, p = 81, B = 10, seed = 123)
}
\arguments{
\item{D}{The dimension of the latent space. Default is 2.}
\item{p}{The number of nodes in the graph. Default is 81.}
\item{B}{The number of conditions. Default is 10.}
\item{seed}{The random seed. Default is 123.}
}
\value{
A list with the following elements:
\item{data}{A list of B elements, where each element contains n-by-p matrix of simulated Gaussian data.}
\item{cloc_true}{A matrix of the true condition-specific intercepts. Dimension is ncol = B, nrow = 2*D.}
\item{alpha_true}{A vector of the true condition-specific slopes. Length is B.}
\item{beta_true}{A scalar of the true slope of the covariate.}
}
\description{
This function simulates data from a latent space model with binary edges and a continuous covariate. The latent space model includes two condition-specific intercepts and a condition-specific slope.
}
\examples{
sim_data <- sim.rgm()
}
