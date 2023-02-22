\name{rgm}
\alias{rgm}
\title{Latent Space Model for Multiple Binary Networks}
\description{
  The function implements the Latent Space Model (LSM) for multiple binary networks. The LSM assumes that there exists an underlying latent space where nodes are located. The model estimates the latent positions, the intercepts for each network, and the network structures for each network.
}
\usage{
  rgm(data, Z=NULL, initial.graphs=NULL, D=2, initial.cloc=NULL, initial.alpha=NULL, initial.beta=NULL, bd.iter=100,iter=1000,burnin=NULL, method = c("ggm","gcgm"), gcgm.dwpar=NULL)
}
\arguments{
  \item{data}{a list of B binary network adjacency matrices, where B is the number of networks}
  \item{Z}{a design matrix for the covariates. Default is \code{NULL}}
  \item{initial.graphs}{an optional list of B binary adjacency matrices for the initial graphs. Default is \code{NULL}}
  \item{D}{number of dimensions in the latent space. Default is 2}
  \item{initial.cloc}{initial values for the B x D matrix of latent node positions. Default is \code{NULL}}
  \item{initial.alpha}{initial values for the B-dimensional intercepts. Default is \code{NULL}}
  \item{initial.beta}{initial values for the regression coefficients for the covariates in Z. Default is \code{NULL}}
  \item{bd.iter}{number of iterations for the BDgraph function. Default is 100}
  \item{iter}{number of iterations for the MCMC sampler. Default is 1000}
  \item{burnin}{number of burn-in iterations to discard. If not specified, the function discards the first 75\% of the iterations. Default is \code{NULL}}
  \item{method}{method used for network estimation. Options are "ggm" (Gaussian graphical model) or "gcgm" (Gaussian copula graphical model). Default is "ggm"}
  \item{gcgm.dwpar}{a list of B elements, each containing the parameters for the truncated data in the GC-GM. If method is "gcgm", this parameter is required. Default is \code{NULL}}
}
\details{
The LSM for multiple binary networks is a Bayesian model that infers the location of each node in a D-dimensional latent space and the probability of a link between two nodes as a function of their Euclidean distance in the latent space. The model also allows for network-specific intercepts and regression coefficients for covariates in Z.

The function first initializes the latent positions and intercepts, and the initial graph if not provided. It then loops through the iterations and updates the latent positions and intercepts using the Gmcmc function. Next, it updates the probability of edge connections for each condition and each edge. Finally, it updates the network structure for each condition using the BDgraph package, if method is "ggm". If method is "gcgm", the function updates the data using the sample.data function and tpoints from the truncated data. The function returns the posterior samples of the parameters after discarding the burn-in period.
}
\value{
  A list containing the posterior samples of the model parameters. The list includes:
  \item{sample.alpha}{a B x (iter - burnin) matrix of the \code{sample.alpha} posterior samples of the network-specific intercepts}
  \item{sample.beta}{a matrix with ncol(Z) rows and (iter - burnin) columns of the posterior samples of the regression coefficients for the covariates in Z. This is only returned if Z is not \code{NULL}.}
  \item{sample.cloc}{a B x D x (iter - burnin) array of the posterior samples of the latent positions of the nodes.}
  \item{sample.graphs}{an n.edge x B x (iter - burnin) array of the posterior samples of the network structures.}
  \item{pi.edgpost}{an n.edge x B matrix of the posterior edge probabilities in each network.}
  \item{pi.probit}{an n.edge x B x (iter - burnin) array of the estimated probit probabilities of the edge connections in each network.}
}
\examples{
set.seed(123)
n <- 50
p <- 10
B <- 3
data <- lapply(1:B, function(i) matrix(rbinom(np, 1, 0.2), n, p))
Z <- matrix(rnorm(nB), nrow=n, ncol=B)
run the model

results <- lsm.bd(data = data, Z = Z, iter = 500, bd.iter = 20)
extract the posterior samples

sample.alpha <- results$sample.alpha
sample.beta <- results$sample.beta
sample.cloc <- results$sample.cloc
sample.graphs <- results$sample.graphs
pi.edgpost <- results$pi.edgpost
pi.probit <- results$pi.probit
}
\author{
Veronica Vinciotti, Ernst C. Wit and Francisco Richter
}
