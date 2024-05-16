library(VSURF)
library(Boruta)
library(varSelRF)
library(vita)
library(CoVVSURF)
library(progress)
library(randomForest)
library(mlbench)

simulation <- function(N = 100, P = 5000, p = 6, M = c(10,50), R = 100){
  
  simulations = list()
  k = 1
  
  # Scénarios dépendants
  for (m in M){
    data = list()
    
    for (r in 1:R){
      X <- matrix(runif(p*N), ncol = p)
      epsilon <- rnorm(N, 0, 0.2) 
      Y <- 0.25*exp(4*X[,1]) + 4/(1+exp(-20*(X[,2]-0.5))) + 3*X[,3] + epsilon
      
      for (i in 1:p) {
        v <- matrix(NA, N, m)
        for (j in 1:m){
          delta <- rnorm(N, 0, 0.3)
          v[,j] <- X[,i] + (0.01 + (0.5 * (j - 1)) / (m - 1))*delta
        } 
        X <- cbind(X, v)
      }
      
      var_sup <- matrix(runif(N*(P-p-p*m)), ncol = P-p-p*m)
      X <- cbind(X, var_sup) 
      data[[r]] <- list(x = X, y = Y)
    }
    
    simulations[[k]] <- data
    k <- k+1
  }
  
  # Scénarios indépendants (ou nuls)
  # for (m in M){
  #   data = list()
  # 
  #   for (r in 1:R){
  #     X <- matrix(runif(p*N), ncol = p)
  #     X_ind <- matrix(rnorm(3*N, 0, 0.2), ncol = 3)
  #     epsilon <- rnorm(N, 0, 0.2)
  #     Y <- 0.25*exp(4*X_ind[,1]) + 4/(1+exp(-20*(X_ind[,2]-0.5))) + 3*X_ind[,3] + epsilon
  # 
  #     for (i in 1:p) {
  #       v <- matrix(NA, N, m)
  #       for (j in 1:m){
  #         delta <- rnorm(N, 0, 0.3)
  #         v[,j] <- X[,i] + (0.01 + (0.5 * (j - 1)) / (m - 1))*delta
  #       }
  #       X <- cbind(X, v)
  #     }
  # 
  #     var_sup <- matrix(runif(N*(P-p-p*m)), ncol = P-p-p*m)
  #     X <- cbind(X, var_sup)
  #     data[[r]] <- list(x = X, y = Y)
  #   }
  # 
  #   simulations[[k]] <- data
  #   k <- k+1
  # }
  
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
    
    # vsurf_vsi <- list()
    # vsurf_vsp <- list()
    boruta_vs <- list()
    janitza_vs1 <- list()
    janitza_vs2 <- list()
    # altmann_vs1 <- list()
    # altmann_vs2 <- list()
    
    i = 1
    
    for (replica in simulation){
      
      # VSURF
      # vsurf <- VSURF(replica$x, replica$y)
      # vsurf_vsi[[i]] <- vsurf$varselect.interp
      # vsurf_vsp[[i]] <- vsurf$varselect.predict
      
      # Boruta 
      boruta <- Boruta(replica$x, replica$y)
      boruta_vs[[i]] <- which(boruta$finalDecision=='Confirmed')

      # Janitza
      PerVarImp1 <- CVPVI(replica$x, replica$y)
      janitza <- NTA(PerVarImp1$cv_varim)
      janitza_vs1[[i]] <- which(janitza$pvalue==0)
      janitza_vs2[[i]] <- which(janitza$pvalue<median(janitza$pvalue))
      
      # Altmann
      # RF <- randomForest(replica$x, replica$y, importance = T)
      # PerVarImp2 <- PIMP(replica$x, replica$y, rForest = RF, S = 20)
      # altmann <- PimpTest(PerVarImp2)
      # altmann_vs1[[i]] <- which(altmann$pvalue<0.001)
      # altmann_vs2[[i]] <- which(altmann$pvalue<floor(median(altmann$pvalue)*1000)*0.001)
      
      # CoV/VSURF
      # covsurf <- covsurf(replica$x, replica$y)
      # covsurf$vsurf_ptree$varselect.interp 
      
      
      # Calcul de performances 
      # ICI
      
      i <- i+1
      pb$tick()
    }
    
    # for (replica in simulation[51:100]){
    #   
    # }
    
    j = 1
    vars_select <- list()
    methods <- list(
                    # vsurf_vsi, 
                    # vsurf_vsp, 
                    boruta_vs, 
                    janitza_vs1,
                    janitza_vs2 
                    # altmann_vs1, 
                    # altmann_vs2
                    )
    
    for (method in methods){
      vars <- numeric(5000)
      for (vs in method) { vars[vs] <- vars[vs] + 1 }
      vars_select[[j]] <- which(vars/length(method) > 0.95) 
      j <- j+1
    }
    
    resultats[[k]] <- list(
                           # vsurf_vsi = vars_select[[1]], 
                           # vsurf_vsp = vars_select[[2]], 
                           boruta_vs = vars_select[[1]], 
                           janitza_vs1 = vars_select[[2]],
                           janitza_vs2 = vars_select[[3]] 
                           # altmann_vs1 = vars_select[[6]], 
                           # altmann_vs2 = vars_select[[7]]
                           )
    k <- k+1
    pb$tick()
  }
  
  return(resultats)
}

evaluation2 <- function(simulations){
  
  pb <- progress_bar$new(
    format = "[:bar] :percent ETA: :eta",
    total = length(simulations)*(length(simulations[[1]])/2)+1
  )
  
  # Critères d'évaluation
  sensibilite_tot <- list()
  fdr_tot <- list()
  stabilite_tot <- list()
  resultats <- list()
  
  k = 1
  pb$tick()
  
  # Pour chaque scénario...
  for (simulation in simulations){
    
    methodes <- list(boruta_vs = c(),
                     # vsurf_vsi = c(),
                     # vsurf_vsp = c(),
                     # altmann_vs1 = c(),
                     # altmann_vs2 = c()),
                     janitza_vs = c(),
                     janitza_vs2 = c())
                      
    # Initialisation des critères d'évaluation
    sensibilite <- matrix(nrow = length(simulation)/2, ncol = length(methodes))
    fdr <- matrix(nrow = length(simulation)/2, ncol = length(methodes))
    stabilite <- matrix(nrow = length(simulation)/2-1, ncol = length(methodes))
    
    i = 1
    
    # Pour chaque réplica...
    for (replica in simulation[1:(length(simulation)/2)]){

      # VSURF
      # vsurf <- VSURF(replica$x, replica$y)
      # methodes$vsurf_vsi <- c(vsurf_vsi, vsurf$varselect.interp)
      # methodes$vsurf_vsp <- c(vsurf_vsp, vsurf$varselect.predict)
      
      # Boruta 
      boruta <- Boruta(replica$x, replica$y)
      methodes$boruta_vs <- which(boruta$finalDecision=='Confirmed')

      # Janitza
      PerVarImp1 <- CVPVI(replica$x, replica$y)
      janitza <- NTA(PerVarImp1$cv_varim)
      methodes$janitza_vs <- which(janitza$pvalue==0)
      methodes$janitza_vs2 <- which(janitza$pvalue<0.001)

      # Altmann
      # RF <- randomForest(replica$x, replica$y, importance = T)
      # PerVarImp2 <- PIMP(replica$x, replica$y, rForest = RF, S = 20)
      # altmann <- PimpTest(PerVarImp2)
      # methodes$altmann_vs1 <- which(methodes$altmann$pvalue<0.001)
      # methodes$altmann_vs2 <- which(methodes$altmann$pvalue<floor(median(altmann$pvalue)*1000)*0.001)
      
      # CoV/VSURF
      # covsurf <- covsurf(replica$x, replica$y)
      # covsurf$vsurf_ptree$varselect.interp 
      
      # Sensibilité & FDR
      j <- 1
      for (vars in methodes){
        if (k == 1){
          p <- c(c(1,2,3), seq(7,36,1))
        } else {
          p <- c(c(1,2,3), seq(7,156,1))
        }
        vp <- length(intersect(vars, p))
        fn <- length(p) - vp
        fp <- length(setdiff(vars, p))
        
        sensibilite[i,j] <- vp/(vp+fn)
        fdr[i,j] <- fp/(fp+vp)
        j <- j+1
      }
      
      # Stabilité
      if (i!=1){
        for (v in seq_along(methodes)){
          stabilite[i-1,v] <- length(intersect(methodes[[v]], ancienne_methodes[[v]])) / length(union(methodes[[v]], ancienne_methodes[[v]]))
        }
      }
      ancienne_methodes <- methodes
      
      i <- i+1
      pb$tick()
    }
    
    # Moyenne des critères pour l'ensemble des réplicas d'un scénario
    sensibilite_tot[[k]] <- colMeans(sensibilite)
    fdr_tot[[k]] <- colMeans(fdr)
    stabilite_tot[[k]] <- colMeans(stabilite)
    k <- k+1
    
    # for (replica in simulation[(length(simulation)/2+1):length(simulation)]){
    #   
    # }
  }
  
  resultats <- list(sensibilité = sensibilite_tot, fdr = fdr_tot, stabilité = stabilite_tot)
  return(resultats)
}


simu <- simulation(R=100)
res2 <- evaluation2(simu)
res2$sensibilité
res2$fdr
res2$stabilité




# Test des méthodes individuellement 

## BORUTA ##
test0 <- Boruta(simu[[2]][[2]]$x, simu[[2]][[2]]$y)
which(test0$finalDecision=='Confirmed')

## JANITZA ##
test1 <- CVPVI(simu[[2]][[2]]$x, simu[[2]][[2]]$y) 
test2 <- NTA(test1$cv_varim)
# Visualisation
mean(test2$pvalue)
median(test2$pvalue)
plot(sort(test2$pvalue)[1:5000])
# Valeur seuil
which(test2$pvalue==0)
which(test2$pvalue<median(test2$pvalue))

## ALTMANN ##
RF <- randomForest(simu[[1]][[5]]$x, simu[[1]][[5]]$y, importance = TRUE)
test3 <- PIMP(simu[[1]][[5]]$x, simu[[1]][[5]]$y, RF)
test4 <- PimpTest(test3)
mean(test4$pvalue)
median(test4$pvalue)
plot(sort(test4$pvalue))
which(test4$pvalue==0)

# VSURF 
test5 <- VSURF(simu[[2]][[1]]$x, simu[[2]][[1]]$y)
sort(test5$varselect.interp)
test5$mean.perf

# COVSURF
test6 <- covsurf(simu[[1]][[5]]$x, simu[[1]][[5]]$y, kval = c(2:3))
sort(test6$vsurf_ptree$varselect.interp)

