sim.rgm <- function(n=346, D = 2, p = 87, B = 13, seed = 123, mcmc_iter=100, alpha=NULL,theta=NULL,loc=NULL, X=NULL) {

  set.seed(seed)

  # Simulating latent space model
  if(is.null(loc)){
    cond.loc1 <- rnorm(B, 0, 0.3)
    cond.loc2 <- rnorm(B, 0, 0.3)
    cloc.true <- cbind(cond.loc1, cond.loc2)
  }  else
    cloc.true<-loc

  # Simulating one covariate
  n.edge <- p*(p-1)/2
  if(is.null(X)){
    X <- runif(n.edge, -0.5, 0.5)
    X <- as.matrix(X)
  } else
    X<-as.matrix(X)

  # True parameters
  if(is.null(alpha))
    alpha.true <- rnorm(B, -2)
  else
    alpha.true<-alpha
  if(is.null(theta))
     beta.true<-2.5
   else
    beta.true <- theta

  # Simulating true graph
  G.true <- matrix(0, ncol = B, nrow = n.edge)
  dist.cond <- matrix(ncol=B, nrow=n.edge)
  m <- matrix(1:p, ncol=p, nrow=p)
  e1 <- t(m)[lower.tri(m)]
  e2 <- m[lower.tri(m)]
  Pi.true <- matrix(ncol = B, nrow = n.edge)
  sp<-NULL
  for(i in 1:mcmc_iter) {
    for (b in 1:B){
      # Updating condition-specific intercept
     dist.cond[,b] <- apply(G.true, 1, function(g, cloc, b) { crossprod(apply(cloc*g, 2, sum) - cloc[b,]*g[b], cloc[b,]) }, cloc = cloc.true, b = b)
     for (k in 2:p){
      for (j in 1:(k-1)){
        ind<-e1==j & e2==k
        Pi.true[ind,b]<-pnorm(alpha.true[b]+dist.cond[ind,b]+X[ind,]%*%beta.true)
      }
    }
    G.true[,b] = rbinom(n.edge,1,Pi.true[,b])
  }
  sp<-c(sp,sum(G.true))
}

  # Simulating data (Gaussian)
  data <- vector(mode = "list", length = B)
  for (j in 1:B) {
    A <- matrix(0, nrow = p, ncol = p)
    A[lower.tri(A)] <- G.true[,j]
    A <- A + t(A)
    data[[j]] <- BDgraph::bdgraph.sim( p = p, n = n, graph = A)$data
  }

  list(data = data, X=X, loc = cloc.true, alpha = alpha.true, theta = beta.true,G=G.true,diagnostic=sp)
}
