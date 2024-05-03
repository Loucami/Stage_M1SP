

simulation1 <- function(N = 100, P = 5000, p = 6, M = c(10,50)){
  
  simulations = list()
  k = 1
  
  # Scénarios dépendants
  for (m in M){
    data = list()
    for (r in 1:R){
      X <- matrix(runif(p*N), ncol = p)
      epsilon <- rnorm(N, 0, 0.2)
      Y <- 0.25*exp(4*X[,1]) + 4/(1+exp(-20*(X[,2]-0.5))) + X[,3] + epsilon
      
      for (i in 1:p) {
        v <- matrix(NA,N,m)
        for (j in 1:m){v[,j] <- X[,i] + (0.01 + (0.5 * (j - 1)) / (m - 1))}
        X <- cbind(X, v)
      }
      
      X <- cbind(X,matrix(0,N,p_tot-p-p*m))
      data[[r]] <- list(x = X, y = Y)
    }
    
    simulations[[k]] <- data
    k <- k+1
  }
  
  # Scénarios indépendants (ou nuls)
  for (m in M){
    data = list()
    for (r in 1:R){
      X <- matrix(runif(p*N), ncol = p)
      X_ind <- matrix(rnorm(3*N,0,0.2), ncol = 3)
      epsilon <- rnorm(N, 0, 0.2)
      Y <- 0.25*exp(4*X_ind[,1]) + 4/(1+exp(-20*(X_ind[,2]-0.5))) + X_ind[,3] + epsilon
      
      for (i in 1:p) {
        v <- matrix(NA,N,m)
        for (j in 1:m){v[,j] <- X[,i] + (0.01 + (0.5 * (j - 1)) / (m - 1))}
        X <- cbind(X, v)
      }
      
      X <- cbind(X,matrix(0,N,p_tot-p-p*m))
      data[[r]] <- list(x = X, y = Y)
    }
    
    simulations[[k]] <- data
    k <- k+1
  }
  
  return(simulations)
}


