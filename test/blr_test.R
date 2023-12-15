# generate data
set.seed(1)
n <- 100
X <- matrix(rnorm(n * 2), ncol = 2)
theta_true <- c(1, 2)
z <- X
y <- ifelse(z > 0, 1, 0)
y <- y[,1]
# fit the model

time0 = proc.time()
theta_hat <- blr(y, X, theta = c(0,0), theta_0 = c(0,0), N_sim = 10000,offset = 0)
time1 = proc.time()
time1 - time0

# summary of posterior samples
mean(theta_hat[,1])
mean(theta_hat[,2])


time0 = proc.time()
theta_hat <- blr_fast(y, X, theta = c(0,0), theta_0 = c(0,0), N_sim = 10000,offset = 0)
time1 = proc.time()
time1 - time0

# summary of posterior samples
mean(theta_hat$theta[,1])
mean(theta_hat$theta[,2])

set.seed(123)
n <- 100
p <- 5
X <- matrix(rnorm(n * p), n, p)
theta <- c(2, 1, -1, 0, 0)
offset <- rep(0, n)
mu <- X %*% theta + offset
y <- rbinom(n, 1, plogis(mu))
theta_chain <- blr2(y, X, offset = offset, theta_0 = rep(0, p), N_sim = 10000000)

# summary of posterior samples
mean(theta_chain[,1])
mean(theta_chain[,2])
mean(theta_chain[,3])
mean(theta_chain[,4])
mean(theta_chain[,5])

