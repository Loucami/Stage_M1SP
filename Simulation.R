library(VSURF)
library(Boruta)
library(varSelRF)
library(vita)
library(CoVVSURF)
library(progress)


simulation1 <- function(N = 100, P = 5000, p = 6, M = c(10,50), R = 100){
  
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
      
      X <- cbind(X,matrix(0,N,P-p-p*m))
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
      
      X <- cbind(X,matrix(0,N,P-p-p*m))
      data[[r]] <- list(x = X, y = Y)
    }
    
    simulations[[k]] <- data
    k <- k+1
  }
  
  return(simulations)
}

 
evaluation1 <- function(simulations){
  
  pb <- progress_bar$new(
    format = "[:bar] :percent ETA: :eta",
    total = length(simulations)*(length(simulations[[1]])+1)+1
  )
  
  pb$tick()
  resultats <- list()
  k = 1
  
  for (simulation in simulations){
    
    vsurf_vsi <- list()
    vsurf_vsp <- list()
    boruta_vs <- list()
    
    i = 1
    
    for (replica in simulation){
      
      # VSURF
      # vsurf <- VSURF(replica$x, replica$y)
      # vsurf_vsi[[i]] <- vsurf$varselect.interp
      # vsurf_vsp[[i]] <- vsurf$varselect.predict
      
      # Boruta 
      boruta <- Boruta(replica$x, replica$y)
      boruta_vs[[i]] <- which(boruta$finalDecision=='Confirmed')
      # plot(boruta)

      # Janitza
      # PerVarImp1 <- CVPVI(replica$x, replica$y)
      # janitza <- NTA(PerVarImp1$cv_varim)
      # print(summary(janitza))
      # janitza_vs[[i]] <- which(janitza$pvalue<0.01)
      
      # Altmann
      # PerVarImp2 <- PIMP(replica$x, replica$y, randomForest(replica$x, replica$y))
      # altmann_vs[[i]] <- PimpTest(PerVarImp2)
      # altmann_vs[[i]] <- which(altmann_vs$pvalue<0.05)
      
      # CoV/VSURF
      # covsurf <- covsurf(replica$x, replica$y)
      # covsurf$vsurf_ptree$varselect.interp 
      
      i <- i+1
      pb$tick()
    }
    
    # for (replica in simulation[51:100]){
    #   
    # }
    
    j = 1
    vars_select <- list()
    methods <- list(boruta_vs)
    for (method in methods){
      vars <- numeric(5000)
      for (vs in method) { vars[vs] <- vars[vs] + 1 }
      vars_select[[j]] <- which(vars/length(method) > 0.95) 
      j <- j+1
    }
    
    resultats[[k]] <- list(boruta = vars_select[[1]])
    k <- k+1
    pb$tick()
  }
  
  return(resultats)
}

simu <- simulation1(R=1)
res <- evaluation1(simu)
