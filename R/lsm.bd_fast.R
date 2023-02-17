lsm.bd_fast<-function(data,X=NULL,iter=1000,burnin=0,initial.graphs=NULL, D=2, initial.loc=NULL, initial.alpha=NULL, initial.theta=NULL, bd.iter=20,bd.jump=10, method=c("ggm","gcgm"), gcgm.dwpar=NULL)
{
  p<-ncol(data[[1]]) #number of nodes
  n.edge<-p*(p-1)/2 #number of edges
  B<-length(data) #number of conditions
  
  sample.graphs<-array(dim=c(n.edge,B,iter))
  sample.K<-array(dim=c(n.edge,B,iter))
  sample.pi<-array(dim=c(n.edge,B,iter))
  
  # lower triangle of pxp matrix
  m<-matrix(1:p,ncol=p,nrow = p)
  lt<-lower.tri(m)

  if(is.null(initial.graphs))
  {
    for(i in 1:B)
    {
      g<-huge.select(huge(as.matrix(data[[i]]),method="glasso"),criterion="stars")$refit
      sample.graphs[,i,1]<-g[lt]
    }
  }
  else
    sample.graphs[,,1]<-initial.graphs
  
  sample.cloc<-array(dim = c(B,D,iter))
  sample.alpha<-matrix(0,B,iter)
  
  if(is.null(initial.loc))
    sample.cloc[,,1]<-matrix(rnorm(B*D),ncol=D)
  else
    sample.cloc[,,1]<-initial.loc
  
  Z<-X
  if(!is.null(X))
  {
    Z<-as.matrix(X)
    sample.beta<-matrix(0,ncol(Z),iter)
    if(is.null(initial.theta))
    {
      y<-as.vector(sample.graphs[,,1])
      X<-apply(Z,2,rep,B)
      sample.beta[,1]<-coef(glm(y~X, family=binomial(link = "probit")))[-1]
    }
    if(!is.null(initial.theta))
      sample.beta[,1]<-initial.theta
  }
  
  if(is.null(burnin))
    burnin<-floor(0.75*iter)
  
  # edge indicators
  e1<-t(m)[lt]
  e2<-m[lt]
  
    
  
  # Initialize K (precision matrix)
  K=vector(mode="list", length=B)
  for (i in 1:B){
    K[[i]]<-diag(p)
    sample.K[,i,1]<-K[[i]][lt]
  }
  
  
  if (method[1]=="gcgm"){
    discrete.data<-data
    #calculate truncated points
    tpoints<-vector("list",B)
    for(i in 1:B)
    {
      tpoints[[i]]<-vector("list",2)
      beta.dw<-gcgm.dwpar[[i]]$beta
      q<-gcgm.dwpar[[i]]$q
      pii<-matrix(rep(gcgm.dwpar[[i]]$pii,each=nrow(q)),nrow(q),ncol(q))
      pdw_lb = BDgraph::pdweibull( data[[i]] - 1, q = q, beta = beta.dw)
      pdw_ub = BDgraph::pdweibull( data [[i]], q = q, beta = beta.dw)
      tpoints[[i]][[1]]<-stats::qnorm( ( 1 - pii)*( data[[i]] != 0 ) + pii*pdw_lb)
      tpoints[[i]][[2]] <- stats::qnorm( (1 - pii) + pii*pdw_ub)
    }
  }
  
  pb <- txtProgressBar(min = 0, max = (iter-1), style = 3)
  for (k in 1: (iter-1))
  {
    setTxtProgressBar(pb = pb, value = k,title = "Performing MCMC iterations") 

    # update data if the Gaussian Copula GM (gcgm) is selected
    if (method[1]=="gcgm"){
      data<-sample.data(data,discrete.data, K, tpoints)
    }
    
    
    # update latent node and condition locations
    G<-sample.graphs[,,k]
    if(is.null(Z))
      G.loc<-Gmcmc_fast(G,alpha=sample.alpha[,k],cloc=sample.cloc[,,k],n.iter=1,n.burnin = 0)
    else
      G.loc<-Gmcmc_fast(G,Z=Z,alpha=sample.alpha[,k],beta=sample.beta[,k],cloc=sample.cloc[,,k],n.iter=1,n.burnin = 0)
    
    cloc<-G.loc$cloc[,,1]
    alpha<-G.loc$alpha
    beta<- G.loc$beta
    
    dist.cond<-matrix(ncol=B,nrow=n.edge)
    for (b in 1:B){
      #updating condition-specific intercept
      dist.cond[,b]<-apply(G,1,function(g,cloc,b){crossprod(colSums(cloc * g)-cloc[b,]*g[b],cloc[b,])},cloc=cloc,b=b)
    }
    Pi = matrix(ncol=B,nrow=n.edge)
    for (b in 1:B){
      for (i in 2:p){
        for (j in 1:(i-1)){
          ind<-e1==j & e2==i
          if(is.null(Z))
            Pi[ind,b]<-pnorm(alpha[b]+dist.cond[ind,b])
          else
            Pi[ind,b]<-pnorm(alpha[b]+dist.cond[ind,b]+Z[ind,]%*%beta)
        }
      }
    }

    for (j in 1:B)
    {
      pi.post<-Pi[,j]
      g.prior<-matrix(0,nrow=p,ncol=p)
      g.prior[lt] <- pi.post
      g.prior<-g.prior+t(g.prior)
      
      g.start<-matrix(0,nrow=p,ncol=p)
      g.start[lt] <- sample.graphs[,j,k]
      g.start<-g.start+t(g.start)
      
      # update K
      res.bd<-BDgraph::bdgraph(data[[j]], iter = bd.iter, jump=bd.jump, g.start=g.start,  g.prior=g.prior, save=FALSE, burnin=0,verbose = FALSE)
     
      g<-res.bd$last_graph
      K[[j]]<-res.bd$last_K
      sample.graphs[,j,k+1]<-g[lt]
      sample.K[,j,k+1] <- K[[j]][lt]
      sample.pi[,j,k+1] <- t(res.bd$p_links)[lt]
    }
    
    sample.cloc[,,k+1]<-cloc
    sample.alpha[,k+1]<-alpha
    if(!is.null(Z))
      sample.beta[,k+1]<-beta
  }
  
  sample.cloc<-sample.cloc[,,-(1:burnin)]
  sample.alpha<-sample.alpha[,-(1:burnin)]
  sample.graphs<-sample.graphs[,,-(1:burnin)]
  sample.K<-sample.K[,,-(1:burnin)]
  sample.pi<-sample.pi[,,-(1:burnin)]
  
  if(!is.null(Z))
    sample.beta<-sample.beta[,-(1:burnin),drop=FALSE]
  
  ##probit probabilities from latent space
  n.iter<-dim(sample.cloc)[3]
  pi.probit = array(dim=c(n.edge,B,n.iter))
  for (k in 1: n.iter){
    G<-sample.graphs[,,k]
    alpha<-sample.alpha[,k]
    if(!is.null(Z))
      beta<-sample.beta[,k]
    cloc<-sample.cloc[,,k]
    dist.cond<-matrix(ncol=B,nrow=n.edge)
    for (b in 1:B){
      #updating condition-specific intercept
      dist.cond[,b]<-apply(G,1,function(g,cloc,b){crossprod(apply(cloc*g,2,sum)-cloc[b,]*g[b],cloc[b,])},cloc=cloc,b=b)
    }
    Pi = matrix(ncol=B,nrow=n.edge)
    for (b in 1:B){
      for (i in 2:p){
        for (j in 1:(i-1)){
          ind<-e1==j & e2==i
          if(is.null(Z))
            pi.probit[ind,b,k]<-pnorm(alpha[b]+dist.cond[ind,b])
          if(!is.null(Z))
            pi.probit[ind,b,k]<-pnorm(alpha[b]+dist.cond[ind,b]+Z[ind,]%*%beta)
        }
      }
    }
  }
  if(is.null(Z))
    return(list(sample.alpha=sample.alpha,sample.loc=sample.cloc,sample.K=sample.K,sample.graphs=sample.graphs,sample.pi=sample.pi,pi.probit=pi.probit))
  else
    return(list(sample.alpha=sample.alpha,sample.theta=sample.beta,sample.K=sample.K,sample.loc=sample.cloc,sample.graphs=sample.graphs,sample.pi=sample.pi,pi.probit=pi.probit))
}


