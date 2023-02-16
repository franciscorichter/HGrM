lsm.bd_fast<-function(data,Z=NULL,initial.graphs=NULL, D=2, initial.cloc=NULL, initial.alpha=NULL, initial.beta=NULL, bd.iter=100,iter=1000,burnin=NULL, gcgm.dwpar=NULL)
{
  p<-ncol(data[[1]]) #number of nodes
  n.edge<-p*(p-1)/2 #number of edges
  B<-length(data) #number of conditions
  pi.edgpost<-matrix(0,n.edge,B)
  khat.edgpost=NULL
  
  sample.graphs<-array(dim=c(n.edge,B,iter))
  if(is.null(initial.graphs))
  {
    for(i in 1:B)
    {
      g<-huge.select(huge(as.matrix(data[[i]]),method="glasso"),criterion="stars")$refit
      sample.graphs[,i,1]<-g[lower.tri(g)]
    }
  }
  else
    sample.graphs[,,1]<-initial.graphs
  
  sample.cloc<-array(dim = c(B,D,iter))
  sample.alpha<-matrix(0,B,iter)
  
  if(is.null(initial.cloc))
    sample.cloc[,,1]<-matrix(rnorm(B*D),ncol=D)
  else
    sample.cloc[,,1]<-initial.cloc
  
  if(!is.null(Z))
  {
    Z<-as.matrix(Z)
    sample.beta<-matrix(0,ncol(Z),iter)
    if(is.null(initial.beta))
    {
      y<-as.vector(sample.graphs[,,1])
      X<-apply(Z,2,rep,B)
      sample.beta[,1]<-coef(glm(y~X, family=binomial(link = "probit")))[-1]
    }
    if(!is.null(initial.beta))
      sample.beta[,1]<-initial.beta
  }
  log.graph.prob<-NULL
  
  if(is.null(burnin))
    burnin<-floor(0.75*iter)
  
  m<-matrix(1:p,ncol=p,nrow = p)
  e1<-t(m)[lower.tri(m)]
  e2<-m[lower.tri(m)]
  
  
  # Initialize K (precision matrix)
  K=vector(mode="list", length=B)
  for (i in 1:B){
    K[[i]]<-diag(p)
  }
  
  
  pb <- txtProgressBar(min = 0, max = (iter-1), style = 3)
  for (k in 1: (iter-1))
  {
    setTxtProgressBar(pb = pb, value = k,title = "Performing MCMC iterations") 

    # update latent node and condition locations
    G<-sample.graphs[,,k]
    if(is.null(Z))
      G.loc<-Gmcmc_fast(G,alpha=sample.alpha[,k],cloc=sample.cloc[,,k],n.iter=1,n.burnin = 0)
    else
      G.loc<-Gmcmc_fast(G,Z=Z,alpha=sample.alpha[,k],beta=sample.beta[,k],cloc=sample.cloc[,,k],n.iter=1,n.burnin = 0)
    
    cloc<-G.loc$cloc[,,1]
    alpha<-G.loc$alpha
    beta<- G.loc$beta
    log.graph.prob<-c(log.graph.prob,G.loc$log.graph.prob)
    
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
    K<-vector("list",B)
    for (j in 1:B)
    {
      pi.post<-Pi[,j]
      g.prior<-matrix(0,nrow=p,ncol=p)
      g.prior[lower.tri(g.prior)] <- pi.post
      g.prior<-g.prior+t(g.prior)
      
      g.start<-matrix(0,nrow=p,ncol=p)
      g.start[lower.tri(g.start)] <- sample.graphs[,j,k]
      g.start<-g.start+t(g.start)
      if(k < (iter-1))
        # update K
        res.bd<-BDgraph::bdgraph(data[[j]], iter = bd.iter, g.start=g.start,  g.prior=g.prior, save=FALSE, burnin=0,verbose = FALSE)
      else
      {
        res.bd<-BDgraph::bdgraph(data[[j]], iter = bd.iter*100, g.start=g.start,  g.prior=g.prior, save=TRUE, burnin=0, verbose = FALSE)
        pp<-plinks(res.bd)
        pi.edgpost[,j]<-t(pp)[lower.tri(pp)]
        khat.edgpost[[j]]<-res.bd$K_hat[lower.tri(res.bd$K_hat)]
        
      }
      g<-res.bd$last_graph
      K[[j]]<-res.bd$last_K
      sample.graphs[,j,k+1]<-g[lower.tri(g)]
      
      
    }
    sample.cloc[,,k+1]<-cloc
    sample.alpha[,k+1]<-alpha
    if(!is.null(Z))
      sample.beta[,k+1]<-beta
    
  }
  
  sample.cloc<-sample.cloc[,,-(1:burnin)]
  sample.alpha<-sample.alpha[,-(1:burnin)]
  sample.graphs<-sample.graphs[,,-(1:burnin)]
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
    return(list(sample.alpha=sample.alpha,sample.cloc=sample.cloc,khat.edgpost=khat.edgpost,sample.graphs=sample.graphs,pi.edgpost=pi.edgpost,pi.probit=pi.probit))
  else
    return(list(sample.alpha=sample.alpha,sample.beta=sample.beta,khat.edgpost=khat.edgpost,sample.cloc=sample.cloc,sample.graphs=sample.graphs,pi.edgpost=pi.edgpost,pi.probit=pi.probit))
}


