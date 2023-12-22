sample.data<-function(data, K, tpoints){
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


post_processing_rgm <- function(simulated_data,results){
  ## Data .  Extracting data from simulation object
  a = simulated_data
  res = results
  alpha.true<-a$alpha
  beta.true<-a$theta
  cloc.true<-a$loc
  G.true<-a$G
  data<-a$data
  X<-a$X
  ## Estimation
  iter<-ncol(res$sample.alpha)
  #Extracting samples after burnin
  burn<-floor(0.75*iter)
  sample.graphs<-res$sample.graphs[,,-(1:burn)]
  sample.cloc<-res$sample.loc[,,-(1:burn)]
  sample.alpha<-res$sample.alpha[,-(1:burn)]
  sample.beta<-res$sample.theta[,-(1:burn)]
  post.pi<-res$sample.pi[,,-(1:burn)]
  probit.pi<-res$pi.probit[,,-(1:burn)]

  #Applying rotation of latent coordinates
  hlp<-array(apply(sample.cloc,3,rot),dim=dim(sample.cloc))
  sample.cloc<-hlp

  #Mean posterior estimates of the parameters
  cloc.est<-apply(sample.cloc,c(1,2),mean)
  alpha.est<-apply(sample.alpha,1,mean)
  sample.beta<-t(as.matrix(sample.beta))
  beta.est<-apply(sample.beta,1,mean)

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
    geom_point(size = 2, shape = 15) +
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
    label = c("True", "Estimation")
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
    labs(title = "Beta Convergence Over Iterations")


  list(rgm_recovery=rgm_recovery,
       estimation_of_alpha = estimation_of_alpha,
       posterior_distribution=posterior_distribution,
       beta_convergence = beta_convergence)
}



