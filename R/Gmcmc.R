Gmcmc<-function(G, Z=NULL, n.iter=1000,alpha=NULL,beta=NULL,cloc=NULL,n.burnin=500)
{
  B<-ncol(G)
  n.edge<-nrow(G)
  p<-(sqrt(1+8*n.edge)+1)/2
  m<-matrix(1:p,ncol=p,nrow = p)
  e1<-t(m)[lower.tri(m)]
  e2<-m[lower.tri(m)]
  if(is.null(cloc))
    cloc<-matrix(rnorm(B*2),ncol=2)
  if(is.null(alpha))
    alpha<-rnorm(B)

  dim.cond<-ncol(cloc)
  cloc.save<-array(dim = c(B,ncol(cloc),n.iter-n.burnin))
  alpha.save<-matrix(0,nrow=B,ncol=n.iter-n.burnin)
  # log.graph.prob.save<-rep(0,n.iter-n.burnin)

  if(!is.null(Z))
  {
    Z<-as.matrix(Z)
    if(is.null(beta))
      beta<-as.matrix(rep(0,ncol(Z)))
    else
      beta<-as.matrix(beta)
    beta.save<-matrix(0,nrow=ncol(Z),ncol=n.iter-n.burnin)
  }

  #################################################################

  for (k in 1:n.iter){
    y<-as.vector(G)

    #####################################################
    for (b in 1:B){
      # update the latent condition locations
      X<-apply(cloc,2,rep,each=n.edge)*rep(G[,b],B)
      X[(b-1)*n.edge+(1:n.edge),]<-matrix(apply(G,1,function(g,cloc,b){colSums(cloc * g)-cloc[b,]*g[b]},cloc=cloc,b=b),nrow=n.edge,ncol=dim.cond,byrow=T)
      hlp2<-NULL
      for (bb in 1:B){
        if (bb==b){
          hlp2<-c(hlp2,rep(0,n.edge))
        } else {
          hlp3<-apply(G,1,function(g,cloc,bb,b){crossprod(colSums(cloc * g)-cloc[b,]*g[b]-cloc[bb,]*g[bb],cloc[bb,])},cloc=cloc,bb=bb,b=b)
          hlp2<-c(hlp2,hlp3)
        }
      }
      offset<-hlp2+rep(alpha,each=n.edge)
      if(!is.null(Z))
        offset<-hlp2+rep(alpha,each=n.edge)+ rep(Z%*%beta,B)
      cloc[b,]<-bpr(y,X,offset,theta = cloc[b,],theta_0 = rep(0,dim.cond))
    }

    #####################################

    dist.cond<-matrix(ncol=B,nrow=n.edge)
    for (b in 1:B){
      #updating condition-specific intercept
      dist.cond[,b]<-apply(G,1,function(g,cloc,b){crossprod(colSums(cloc * g)-cloc[b,]*g[b],cloc[b,])},cloc=cloc,b=b)
      offset<-dist.cond[,b]
      if(!is.null(Z))
        offset<-dist.cond[,b]+ Z%*%beta
      y<-G[,b]
      X<-as.matrix(rep(1,length(y)))
      alpha[b]<-bpr(y,X,offset,theta = alpha[b],theta_0 = 0)
    }


    if(!is.null(Z)){
      y<-as.vector(G)
      X<-apply(Z,2,rep,B)
      offset<-c(dist.cond)+rep(alpha,each=n.edge)
      beta<-bpr(y,X,offset,theta = beta,theta_0 = rep(0,length(beta)))
      beta<-t(beta)
    }

    if (k>n.burnin){
      cloc.save[,,k-n.burnin]<-cloc
      alpha.save[,k-n.burnin]<-alpha
      if(!is.null(Z))
        beta.save[,k-n.burnin]<-beta
    }
  }

  #################################################################


  if(is.null(Z))
    return(list(alpha=alpha.save,cloc=cloc.save))
  else
    #    return(list(alpha=alpha.save,beta=beta.save,cloc=cloc.save,log.graph.prob=log.graph.prob.save))
    return(list(alpha=alpha.save,beta=beta.save,cloc=cloc.save))
}
