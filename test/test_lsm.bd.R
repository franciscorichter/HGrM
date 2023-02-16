library(HGrM)
########
#Analysis
set.seed(123)
iter<-10 #outer iterations for LSM
bd.iter<-10 #inner iterations for bdgraph

da=sim_HGrM(mcmc_iter = 10)
data = da$data

p<-ncol(data[[1]]) #number of nodes
n.edge<-p*(p-1)/2 #number of edges
B<-length(data) #number of conditions

initial.graphs<-matrix(nrow=n.edge,ncol=B)
for(i in 1:B){
    g<-huge.select(huge(as.matrix(data[[i]]),method="glasso"),criterion="stars")$refit
    initial.graphs[,i]<-g[lower.tri(g)]
}

save(initial.graphs,file = "initialgraphs.RData")

time0=proc.time()
res<-lsm.bd(data = data,
            D=2, 
            bd.iter=bd.iter, 
            iter=iter,
            initial.graphs = initial.graphs)
time1=proc.time()
time1-time0

#save.image("LSM-simresults-p100.RData")

time0=proc.time()
res2<-lsm.bd_fast(data = data,
            D=2,
            bd.iter=bd.iter, 
            iter=iter,
            initial.graphs = initial.graphs)
time1=proc.time()
time1-time0

library(profvis)
profvis({
  time0=proc.time()
  res2<-lsm.bd_fast(data = data,
                    D=2,
                    bd.iter=bd.iter, 
                    iter=iter,
                    initial.graphs = initial.graphs)
  time1=proc.time()
  time1-time0
})
  
  
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
