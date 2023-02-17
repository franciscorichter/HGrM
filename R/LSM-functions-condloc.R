


sample.data<-function(data,discrete.data, K, tpoints){
  B<-length(data)
  p<-ncol(data[[1]])
  dat<-vector(mode="list", length=B)
  for (i in 1:B){
    for (j in 1:p){
      S<-solve(K[[i]])
      S22i<-solve(S[-j,-j])
      S12<-S[j,-j]
      S11<-S[j,j]
      mu.j<-t(S12)%*%S22i%*%t(data[[i]][,-j])
      var.j<-S11-t(S12)%*%S22i%*%as.matrix(S12)
      data[[i]][,j]<-rtruncnorm(length(mu.j),a=tpoints[[i]][[1]][,j],b=tpoints[[i]][[2]][,j], mean=mu.j, sd=sqrt(var.j))
    }
  }
  return(data)
}


rot<-function(loc){
  x.mn<-apply(loc,2,mean)
  alpha<-atan(x.mn[2]/x.mn[1])+(x.mn[1]<0)*pi
  phi<-pi/2-alpha
  angles<-apply(loc,1,function(x){atan(x[2]/x[1])+(x[1]<0)*pi})
  r<-apply(loc,1,function(x){sqrt(x[1]^2+x[2]^2)})
  new.loc<-r*cbind(sin(phi+angles),cos(phi+angles))
  return(new.loc)
}
