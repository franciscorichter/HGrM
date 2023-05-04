graph_search <- function(Z,w,lambda) {
  n <- dim(Z)[1]
  p <- dim(Z)[2]

  # Initialize G, K, and mu
  #huge to get G and K initial
  #G <- matrix(0, p, p)
  K <- diag(p)
  G <- matrix(rbinom(100, 1, 0.5), ncol = 10)
  G.upper <- upper.tri(G)


  # logistic regression with G.ini
  logit <- glm(G[G.upper] ~ w, family = binomial(link="logit"))
  lik_G <- logLik(logit)

  # collect data statistics
  mu <- colMeans(Z)
  S <- cov(Z)

  # prior parameters for the G wishart
  V <- diag(p)
  nu <- 1

  # these values are used in the joint optimization of l_G(K) + l_K(Z)
  fake_S<- (V+S*n)/(n+nu)
  fake_n<- (n+nu)

  # Calculate the initial log-likelihood
  #log_likelihood_current <- log_likelihood(K, G, Z, theta)
  lik_K <- likelihood_K(K, V, S)
  penalty <-lambda*sum(G[G.upper]!=0)
  tot_lik <- lik_K+lik_G+penalty

  # Start the search for the optimal graph
  improving <- TRUE
  while (improving) {
    improving <- FALSE
    best_lik <- tot_lik
    best_edge<- NA

    # Try adding and removing edges
    for (i in 1:(p - 1)) {
      for (j in (i + 1):p) {
        # Add or remove edge (i, j) and update K
        if (G[i, j] == 0) {
          G_new <- add_edge(G, i, j)
          #action <- "add"

          rownames(G_new) = colnames(G_new) = paste0("V",1:nrow(G_new))

          fcg <- fitConGraph(G_new,fake_S,fake_n)
          K_new <- solve(fcg$Shat)
          lik_K_new <- likelihood_K(K = K_new,n = fake_n,S = fake_S)
          # logistic regression with G new
          logit <- glm(G_new[G.upper] ~ w, family = binomial(link="logit"))
          lik_G <- loglik(logit)
          penalty_new <- penalty + lambda
          tot_lik_new <- lik_K_new + lik_G + penalty_new
        } else {
          G_new <- remove_edge(G, i, j)
          #action <- "remove"
          delta_lik_K <- delta_loglik_K_edge_remove(K,i,j,n)
          lik_K_new <- lik_K + delta_lik_K
          # logistic regression with G new
          logit <- glm(G_new[G.upper] ~ w, family = binomial(link="logit"))
          lik_G <- logLik(logit)
          penalty_new <- penalty - lambda
          tot_lik_new <- lik_K_new + lik_G + lambda
        }

        if (tot_lik_new > best_lik){
          best_lik <- tot_lik_new
          best_edge <- c(i,j)
        }
      }
    }

    if (best_lik > tot_lik) {
      improving <- TRUE
      penalty <- penalty + 2*(.5-G[i,j])*lambda
      G[i,j] <- 1-G[i,j]
      G[j,i] <- 1-G[j,i]
      rownames(G) = colnames(G) = paste0("V",1:nrow(G))

      fcg <- fitConGraph(G,fake_S,fake_n)
      K <- solve(fcg$Shat)
      lik_K <- likelihood_K(K = K,n = fake_n,S = fake_S)
      #lik_K <- -fcg$dev/2
    }
  }
  theta <- glm(G[G.upper] ~ w, family = binomial(link="logit"))$coef

  # Return the optimal graph and precision matrix
  list(G = G, K = K,theta=theta)
}
