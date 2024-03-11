#' @importFrom truncnorm rtruncnorm
#' @importFrom MASS ginv
#' @importFrom mvtnorm rmvnorm
#' @importFrom ggplot2 ggplot geom_point geom_abline scale_color_manual labs theme_minimal theme geom_line geom_density geom_vline geom_hline geom_tile scale_fill_gradient2 aes
#' @importFrom pROC roc auc
#' @importFrom reshape2 melt
#' @importFrom stats as.dendrogram binomial coef dist glm hclust order.dendrogram pnorm rbinom rnorm runif
#' @importFrom utils setTxtProgressBar txtProgressBar
#' @importFrom huge huge.select
#' @importFrom huge huge
#' @importFrom ggplot2 element_blank element_text theme_bw


sample.data<-function (data, K, tpoints)
{
  B <- length(data)
  p <- ncol(data[[1]])
  for (i in 1:B) {
    for (j in 1:p) {
      S <- solve(K[[i]])
      S22i <- solve(S[-j, -j])
      S12 <- S[j, -j]
      S11 <- S[j, j]
      mu.j <- t(S12) %*% S22i %*% t(data[[i]][, -j])
      var.j <- S11 - t(S12) %*% S22i %*% as.matrix(S12)
      data[[i]][, j] <- rtruncnorm(length(mu.j), a = tpoints[[i]][[1]][,j], b = tpoints[[i]][[2]][, j], mean = mu.j,sd = sqrt(var.j))
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

bpr <- function(y, X, offset = 0, theta, theta_0 = c(0, 0, 0), N_sim = 1) {


  # Dimensions of theta
  D <- ncol(X)

  # Number of observations
  n <- length(y)
  N1 <- sum(y)
  N0 <- n - N1

  # Conjugate prior on the coefficients theta ~ N(theta_0, Q_0)
  Q_0 <- diag(10, D)

  # Initialize parameters
  z <- rep(NA, n)

  # Matrix storing samples of the theta parameter
  theta_chain <- matrix(0, nrow = N_sim, ncol = D)

  # ---------------------------------
  # Gibbs sampling algorithm
  # ---------------------------------

  # Compute posterior variance of theta
  prec_0 <- ginv(Q_0)
  V <- ginv(prec_0 + t(X) %*% X)

  for (t in 1:N_sim) {
    # Update Mean of z
    mu_z <- X %*% theta + offset
    # Draw latent variable z from its full conditional: z | theta, y, X
    if (N0 > 0) {
      z[y == 0] <- rtruncnorm(N0, mean = mu_z[y == 0], sd = 1, a = -Inf, b = 0)
    }
    if (N1 > 0) {
      z[y == 1] <- rtruncnorm(N1, mean = mu_z[y == 1], sd = 1, a = 0, b = Inf)
    }

    # Compute posterior mean of theta
    M <- V %*% (prec_0 %*% theta_0 + t(X) %*% (z - offset))
    # Draw variable theta from its full conditional: theta | z, X
    theta <- rmvnorm(1, M, V)

    # Store the theta draws
    theta_chain[t, ] <- theta
  }

  return(theta_chain)
}




Gmcmc<-function(G, X=NULL, iter=1000,alpha=NULL,theta=NULL,loc=NULL, burnin=0)
{
  B<-ncol(G)
  n.edge<-nrow(G)
  p<-(sqrt(1+8*n.edge)+1)/2
  m<-matrix(1:p,ncol=p,nrow = p)
  e1<-t(m)[lower.tri(m)]
  e2<-m[lower.tri(m)]
  if(is.null(loc))
    cloc<-matrix(rnorm(B*2),ncol=2)
  else
    cloc <- loc
  if(is.null(alpha))
    alpha<-rnorm(B)

  dim.cond<-ncol(cloc)
  cloc.save<-array(dim = c(B,ncol(cloc), iter-burnin))
  alpha.save<-matrix(0,nrow=B,ncol= iter-burnin)

  Z <- X

  if(!is.null(Z))
  {
    Z<-as.matrix(Z)
    if(is.null(theta))
      beta<-as.matrix(rep(0,ncol(Z)))
    else
      beta<-as.matrix(theta)
    beta.save<-matrix(0,nrow=ncol(Z),ncol=iter-burnin)
  }

  #################################################################

  for (k in 1:iter){
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

    if (k>burnin){
      cloc.save[,,k-burnin]<-cloc
      alpha.save[,k-burnin]<-alpha
      if(!is.null(Z))
        beta.save[,k-burnin]<-beta
    }
  }

  #################################################################


  if(is.null(Z))
    return(list(alpha=alpha.save,loc=cloc.save))
  else
    return(list(alpha=alpha.save,theta=beta.save,loc=cloc.save))
}



post_processing_rgm <- function(simulated_data, results) {
  ## Data extraction from simulation object
  a = simulated_data
  res = results
  alpha.true <- a$alpha
  beta.true <- a$theta
  cloc.true <- a$loc
  G.true <- a$G
  data <- a$data
  X <- a$X



  ## Estimation
  iter <- ncol(res$sample.alpha)
  # Extracting samples after burnin
  burn <- floor(0.75 * iter)
  sample.graphs <- res$sample.graphs[, , -(1:burn)]
  sample.cloc <- res$sample.loc[, , -(1:burn)]
  sample.alpha <- res$sample.alpha[, -(1:burn)]
  sample.beta <- res$sample.theta[, -(1:burn)]
  post.pi <- res$sample.pi[, , -(1:burn)]
  probit.pi <- res$pi.probit[, , -(1:burn)]

  # Applying rotation of latent coordinates
  hlp <- array(apply(sample.cloc, 3, rot), dim = dim(sample.cloc))
  sample.cloc <- hlp

  # Mean posterior estimates of the parameters
  cloc.est <- apply(sample.cloc, c(1, 2), mean)
  alpha.est <- apply(sample.alpha, 1, mean)
  sample.beta <- t(as.matrix(sample.beta))
  beta.est <- apply(sample.beta, 1, mean)


  #Calculating the true edge probabilities (associated to the true graph and true parameters)
  Pi.true<-G.true
  dist.cond<-G.true
  B<-ncol(G.true)
  p<-ncol(data[[1]])
  m <- matrix(1:p, ncol = p, nrow = p)
  e1 <- t(m)[lower.tri(m)]
  e2 <- m[lower.tri(m)]
  for(b in 1:B)
  {
    dist.cond[,b]<-apply(G.true,1,function(g,cloc,b){crossprod(apply(cloc*g,2,sum)-cloc[b,]*g[b],cloc[b,])},cloc=cloc.true,b=b)
    for (i in 2:p){
      for (j in 1:(i-1)){
        ind<-e1==j & e2==i
        Pi.true[ind,b]<-pnorm(alpha.true[b]+dist.cond[ind,b]+X[ind,]%*%beta.true)
      }
    }}

  #Comparing the true edge probabilities with those estimated by the latent probit model, for each environment
  Pi.mean<-apply(probit.pi,c(1,2),mean)

  true_prob = est_prob = true_alpha = est_alpha = beta_value = value = label = iteration = specificity = sensitivity = Var1 = Var2 = NULL

  # Creating a data frame for ggplot
  data_to_plot <- data.frame(
    true_prob = c(Pi.true),
    est_prob = c(Pi.mean),
    environment = factor(rep(1:B, each = nrow(Pi.mean)))
  )

  # Define the colors for each environment
  col_vector <- c('#e6194b', '#3cb44b', '#ffe119', '#4363d8', '#f58231', '#911eb4', '#f032e6', '#fabebe', '#008080', '#000000', "#00FFFFFF", "#80FF00FF", "#FFFF00FF")
  col_vector <- col_vector[1:B]

  # Create the ggplot
  rgm_recovery = ggplot(data_to_plot, aes(x = true_prob, y = est_prob, color = environment)) +
    geom_point() +
    geom_abline(intercept = 0, slope = 1, lwd = 3) +
    scale_color_manual(values = col_vector) +
    labs(x = "True Probit", y = "Estimated Probit", title = "Random Graph Model Recovery") +
    theme_minimal() +
    theme(legend.title = element_blank())


  #Comparing the true alphas with the estimated ones (i.e., checking whether the sparsity levels of the environment-specific graphs are correctly estimated)

  environments=1:B
  # Create a data frame for plotting
  data_to_plot <- data.frame(
    true_alpha = alpha.true,
    est_alpha = alpha.est,
    environment = factor(environments)  # Make sure this matches your data
  )

  # Create the ggplot
  estimation_of_alpha = ggplot(data_to_plot, aes(x = true_alpha, y = est_alpha, color = environment)) +
    geom_point(size = 4, shape = 15) +
    geom_abline(intercept = 0, slope = 1, lwd = 2) +
    scale_color_manual(values = col_vector) +
    labs(x = expression(paste("True ", alpha)),
         y = expression(paste("Posterior Mean ", alpha)),
         title = expression(paste("Estimation of ", alpha))) +
    theme_minimal() +
    theme(legend.title = element_blank(),
          plot.title = element_text(size = 14),
          axis.title = element_text(size = 12))

  beta_samples <- sample.beta[1, ]  # Extracting samples for the first beta parameter

  # Create a data frame for the density plot
  data_to_plot <- data.frame(beta_value = beta_samples)

  # Create a small data frame for the true and estimated lines
  lines_data <- data.frame(
    value = c(a$theta, beta.est[1]),
    label = c("True", "Estimated")
  )

  posterior_distribution = ggplot() +
    geom_density(data = data_to_plot, aes(x = beta_value), fill = "blue", alpha = 0.5) +
    geom_vline(data = lines_data, aes(xintercept = value, color = label), size = .5, linetype = "longdash") +
    scale_color_manual(values = c("True" = "darkgreen", "Estimation" = "red")) +
    labs(x = expression(beta),
         title = expression(paste(beta, " Posterior Distribution"))) +
    theme_minimal() +
    theme(plot.title = element_text(size = 14),
          axis.title = element_text(size = 12),
          legend.position = "top",
          legend.title = element_blank())

  df_p = data.frame(beta=res$sample.theta[1,], iteration=1:ncol(res$sample.theta))

  beta_convergence = ggplot(df_p) +
    geom_line(aes(x = iteration, y = beta)) +
    theme_bw() +
    geom_hline(yintercept = a$theta) +
    labs(title = "Trace plot of beta across MCMC Iterations")


  postpi.mean <- apply(post.pi, c(1, 2), mean)
  predictor <- postpi.mean
  response <- G.true

  roc_data <- lapply(1:B, function(j) {
    roc_obj <- pROC::roc(response = response[, j], predictor = predictor[, j], levels = c(0, 1), direction = "<")
    data.frame(
      specificity = 1 - roc_obj$specificities,
      sensitivity = roc_obj$sensitivities,
      environment = j
    )
  })

  roc_data <- do.call(rbind, roc_data)


  roc_plot = ggplot(roc_data, aes(x = specificity, y = sensitivity, color = as.factor(environment))) +
    geom_line(size = 1) +
    scale_color_manual(values = col_vector) +
    labs(title = "Graph Recovery", x = "Specificity", y = "Sensitivity", color = "Environment") +
    theme_minimal() +
    theme(legend.position = "right")


  #AUC values
  a<-NULL
  for(j in 1:B)
  {
    a[j]<-as.numeric(auc(response = response[,j], predictor = predictor[,j], levels = c(0, 1),direction="<"))
  }

  #Heatmap of posterior edge probabilities for each environment (blue 0, red 1)
  colnames(postpi.mean)<-as.character(1:B)
  dat.sm<-as.matrix(postpi.mean)


  # Perform hierarchical clustering on both rows (edges) and columns (environments)
  row_dend <- as.dendrogram(hclust(dist(dat.sm)))
  col_dend <- as.dendrogram(hclust(dist(t(dat.sm))))

  # Rearrange data according to the clustering
  ordered_dat.sm <- dat.sm[order.dendrogram(row_dend), order.dendrogram(col_dend)]

  # Reshape the matrix into a long format
  long_dat <- reshape2::melt(ordered_dat.sm)
  long_dat$Var1 <- as.factor(long_dat$Var1) # Environments
  long_dat$Var2 <- as.factor(long_dat$Var2) # Edges

  # Create the heatmap using ggplot2
  gg_heatmap <- ggplot(long_dat, aes(x = Var1, y = Var2, fill = value)) +
    geom_tile() +
    scale_fill_gradient2(low = "blue", mid = "white", high = "red",
                         midpoint = 0.5, # Adjust based on your data range
                         limits = c(0, 1)) +
    theme_minimal() +
    theme(axis.text.x = element_blank(),  # Remove x-axis text
          axis.text.y = element_text(angle = 45, hjust = 1)) +
    labs(x = "Edge", y = "Environment", fill = "Posterior Probability")


    # Return the list including the ggplot2 heatmap
    list(rgm_recovery = rgm_recovery,
         estimation_of_alpha = estimation_of_alpha,
         posterior_distribution = posterior_distribution,
         beta_convergence = beta_convergence,
         roc_plot = roc_plot,
         edge_prob = gg_heatmap)
}



