# Méthodes
library(randomForest)
library(VSURF)
library(Boruta)
library(varSelRF)
library(vita)
library(CoVVSURF)
# Barre de progression
library(progress)
# Simulations
library(mlbench)
# Graphiques
library(ggplot2)
library(gridExtra)

simulation_regression_old <- function(N = 100, P = 5000, p = 6, M = c(10,50), R = 100){
  
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

evaluation_old <- function(simulations){
  
  pb <- progress_bar$new(
    format = "[:bar] :percent ETA: :eta",
    total = length(simulations)*(length(simulations[[1]])/2)+3
  )
  
  # Critères d'évaluation
  col_names <- c('Boruta', 'Janitza', 'Altmann', 'VSURF')
  sensibilite_tot <- matrix(nrow = 2, ncol = length(col_names), dimnames = list(c(1,2), col_names))
  fdr_tot <- matrix(nrow = 2, ncol = length(col_names), dimnames = list(c(1,2), col_names))
  stabilite_tot <- matrix(nrow = 2, ncol = length(col_names), dimnames = list(c(1,2), col_names))
  rmse_tot <- matrix(nrow = 2, ncol = length(col_names), dimnames = list(c(1,2), col_names))
  empower_tot <- list()
  time_tot <- matrix(nrow = 2, ncol = length(col_names), dimnames = list(c(1,2), col_names))
  
  valeurs <- list(sensibilite = list(), 
                  fdr = list(),
                  stabilite = list(),
                  rmse = list())
  
  resultats <- list()
  
  k = 1
  pb$tick()
  
  # Pour chaque scénario...
  for (simulation in simulations){
    
    methodes <- list(boruta_vs = c(),
                     janitza_vs = c(),
                     altmann_vs = c(),
                     vsurf_vs = c())
    
    modeles <- list(boruta_rf = c(),
                    janitza_rf = c(),
                    altmann_rf = c(),
                    vsurf_rf = c())
    
    # Initialisation du nombre de variables d'intérêts
    if (k==1) {p <- c(c(1,2,3), seq(7,36,1))} 
    else if (k==2) {p <- c(c(1,2,3), seq(7,156,1))}
    
    # Initialisation des critères d'évaluation
    sensibilite <- matrix(nrow = length(simulation)/2, ncol = length(methodes))
    fdr <- matrix(nrow = length(simulation)/2, ncol = length(methodes))
    stabilite <- matrix(nrow = ((length(simulation)/2)*(length(simulation)/2-1))/2, ncol = length(methodes))
    rmse <- matrix(nrow = length(simulation)/2, ncol = length(methodes))
    empower <- matrix(0, nrow = length(simulation)/2, ncol = length(methodes))
    time <- matrix(0, nrow = length(simulation)/2, ncol = length(methodes))
    
    vars_select <- list()
    RF <- list()
    
    i = 1
    
    # Application des méthodes sur la première moitié des réplicas
    for (replica in simulation[1:(length(simulation)/2)]){
      
      # Boruta 
      time[i,1] <- system.time({
        boruta <- Boruta(replica$x, replica$y)
        methodes$boruta_vs <- which(boruta$finalDecision=='Confirmed')
        modeles$boruta_rf <- randomForest(x = replica$x[,methodes$boruta_vs], y = replica$y)})[3]
      
      # Janitza
      time[i,2] <- system.time({
        PerVarImp1 <- CVPVI(replica$x, replica$y)
        janitza <- NTA(PerVarImp1$cv_varim)
        methodes$janitza_vs <- which(janitza$pvalue==0)
        modeles$janitza_rf <- randomForest(x = replica$x[,methodes$janitza_vs], y = replica$y)})[3]
      
      # Altmann
      time[i,3] <- system.time({
        rf <- randomForest(replica$x, replica$y, importance = T)
        PerVarImp2 <- PIMP(replica$x, replica$y, rForest = rf, S = 50)
        altmann <- PimpTest(PerVarImp2)
        methodes$altmann_vs <- which(altmann$pvalue==0)
        modeles$altmann_rf <- randomForest(x = replica$x[,methodes$altmann_vs], y = replica$y)})[3]
      
      # VSURF
      time[i,4] <- system.time({
        vsurf <- VSURF(replica$x, replica$y)
        methodes$vsurf_vs <- vsurf$varselect.interp
        modeles$vsurf_rf <- randomForest(x = replica$x[,methodes$vsurf_vs], y = replica$y)})[3]
      
      # CoV/VSURF
      # covsurf <- covsurf(replica$x, replica$y)
      # covsurf$vsurf_ptree$varselect.interp 
      
      # Sensibilité & FDR
      j <- 1
      for (vars in methodes){
        vp <- length(intersect(vars, p))
        fn <- length(p) - vp
        fp <- length(setdiff(vars, p))
        
        sensibilite[i,j] <- vp/(vp+fn)
        fdr[i,j] <- fp/(fp+vp)
        j <- j+1
      }
      
      vars_select[[i]] <- methodes
      RF[[i]] <- modeles
      i <- i+1
      pb$tick()
    }
    
    # Stabilité
    pair_index <- 1
    for (a in 1:(length(vars_select) - 1)) {
      for (b in (a + 1):length(vars_select)) {
        for (m in 1:length(methodes)) {
          vars_a <- vars_select[[a]][[m]]
          vars_b <- vars_select[[b]][[m]]
          intersection_ab <- length(intersect(vars_a, vars_b))
          union_ab <- length(union(vars_a, vars_b))
          stabilite[pair_index, m] <- intersection_ab / union_ab
        }
        pair_index <- pair_index + 1
      }
    }
    
    # RMSE 
    i = 1
    # Application des modèles séletcionnés sur l'autre moitié des réplicas
    for (replica in simulation[(length(simulation)/2+1):length(simulation)]){
      for (m in 1:length(methodes)){
        prediction <- predict(RF[[i]][[m]], newdata = replica$x[,vars_select[[i]][[m]]])
        rmse[i,m] <- sqrt(mean((replica$y-prediction)^2))
      }
      i <- i+1
    }
    
    # Empirical Power
    i <- 1
    empower <- list(boruta_ep = matrix(0, length(simulation)/2, 5000), 
                    janitza_ep = matrix(0,  length(simulation)/2, 5000), 
                    altmann_ep = matrix(0,  length(simulation)/2, 5000),
                    vsurf_ep = matrix(0,  length(simulation)/2, 5000))
    for (v in 1:length(vars_select)){
      for (m in 1:length(methodes)){
        empower[[m]][i,vars_select[[i]][[m]]] <- empower[[m]][i,vars_select[[i]][[m]]] + 1
      }
      i <- i+1
    }
    
    # Enregistrement des valeurs pour les IC
    valeurs$sensibilite[[k]] <- sensibilite
    valeurs$fdr[[k]] <- fdr
    valeurs$stabilite[[k]] <- stabilite
    valeurs$rmse[[k]] <- rmse
    
    # Moyenne des critères pour un scénario
    sensibilite_tot[k,] <- colMeans(sensibilite)
    fdr_tot[k,] <- colMeans(fdr)
    stabilite_tot[k,] <- colMeans(stabilite)
    rmse_tot[k,] <- colMeans(rmse)
    for (e in 1:length(empower)){empower[[e]] <- colMeans(empower[[e]])}
    empower_tot[[k]] <- empower
    time_tot[k,] <- colMeans(time)
    
    k <- k+1
    pb$tick()
  }
  
  # Calcul des intervalles de confiance
  # 1 - Distribution Normale
  c <- 1
  criteres <- list(sensibilite_tot, 
                   fdr_tot, 
                   stabilite_tot, 
                   rmse_tot)
  ic <- list(ic_sensi = data.frame(), 
             ic_fdr = data.frame(), 
             ic_stabi = data.frame(), 
             ic_rmse = data.frame())
  for (critere in criteres){
    for (i in 1:nrow(critere)){
      for (j in 1:ncol(critere)){
        z_score <- qnorm(0.975)
        n <- length(simulations[[1]])/2
        if (c<4){
          inf <- round(critere[i,j] - z_score * sqrt(critere[i,j] * (1 - critere[i,j]) / n), 3)
          inf <- ifelse(inf<0, 0, inf)
          sup <- round(critere[i,j] + z_score * sqrt(critere[i,j] * (1 - critere[i,j]) / n), 3)
          sup <- ifelse(sup>1,1, sup)
          ic[[c]][i,j] <- paste(inf, sup, sep = ';')
        } else {
          inf <- round(critere[i,j] - z_score * (sd(valeurs[[c]][[i]][,j]) / sqrt(n)), 3)
          sup <- round(critere[i,j] + z_score * (sd(valeurs[[c]][[i]][,j]) / sqrt(n)), 3)
          ic[[c]][i,j] <- paste(inf, sup, sep = ';')
        }
        colnames(ic[[c]])[j] <- col_names[j]
      }
    }
    c <- c+1
  }
  
  
  resultats <- list(sensibilité = list(valeurs = sensibilite_tot, ic = ic$ic_sensi), 
                    fdr = list(valeurs = fdr_tot, ic = ic$ic_fdr), 
                    stabilité = list(valeurs = stabilite_tot, ic = ic$ic_stabi), 
                    rmse = list(valeurs = rmse_tot, ic = ic$ic_rmse), 
                    empower = empower_tot,
                    time = time_tot)
  
  return(resultats)
}

simulation_regression <- function(N = 100, P = 5000, p = 6, M = c(10,50), R = 100){
  
  simulations = list()
  k = 1
  
  # Scénarios dépendants
  for (m in M){
    data = list()
    
    for (r in 1:R){
      # On crée p variables (x) afin de générer Y avec x1, x2 et x3
      x <- matrix(runif(p*N), ncol = p)
      epsilon <- rnorm(N, 0, 0.2) 
      Y <- 0.25*exp(4*x[,1]) + 4/(1+exp(-20*(x[,2]-0.5))) + 3*x[,3] + epsilon
      
      # On crée p variables (x) pour générer p*m variables corrélées
      vars <- matrix(NA, N, p*m)
      for (i in 1:p) {
        for (j in 1:m){
          delta <- rnorm(N, 0, 0.3)
          vars[,(i-1)*m+j] <- x[,i] + (0.01 + (0.5 * (j - 1)) / (m - 1))*delta
        } 
      }
      
      # On crée notre échantillon X en agrégeant les variables corrélées (vars) avec des variables indépendantes (vars_sup)
      vars_sup <- matrix(runif(N*(P-p*m)), ncol = P-p*m)
      X <- cbind(vars, vars_sup) 
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
  #     x_ind <- matrix(rnorm(3*N, 0, 0.2), ncol = 3)
  #     epsilon <- rnorm(N, 0, 0.2)
  #     Y <- 0.25*exp(4*x_ind[,1]) + 4/(1+exp(-20*(x_ind[,2]-0.5))) + 3*x_ind[,3] + epsilon
  # 
  #     x <- matrix(runif(p*N), ncol = p)
  #     vars <- matrix(NA, N, p*m)
  #     for (i in 1:p) {
  #       for (j in 1:m){
  #         delta <- rnorm(N, 0, 0.3)
  #         vars[,(i-1)*m+j] <- x[,i] + (0.01 + (0.5 * (j - 1)) / (m - 1))*delta
  #       }
  #     }
  # 
  #     vars_sup <- matrix(runif(N*(P-p*m)), ncol = P-p*m)
  #     X <- cbind(vars, vars_sup) 
  #     data[[r]] <- list(x = X, y = Y)
  #   }
  # 
  #   simulations[[k]] <- data
  #   k <- k+1
  # }
  
  return(simulations)
}

evaluation <- function(simulations){
  
  pb <- progress_bar$new(
    format = "[:bar] :percent ETA: :eta",
    total = length(simulations)*(length(simulations[[1]])/2)+3
  )
  
  # Critères d'évaluation
  col_names <- c('Boruta', 'Janitza')
  empower_tot <- list()
  time_tot <- matrix(nrow = 2, ncol = length(col_names), dimnames = list(c(1,2), col_names))
  valeurs <- list(sensibilite = list(), 
                  fdr = list(),
                  stabilite = list(),
                  rmse = list())
  
  resultats <- list()
  
  k = 1
  pb$tick()
  
  # Pour chaque scénario...
  for (simulation in simulations){
    
    methodes <- list(boruta_vs = c(),
                     janitza_vs = c())
    
    modeles <- list(boruta_rf = c(),
                    janitza_rf = c())
    
    # Initialisation du nombre de variables d'intérêts
    if (k==1) {p <- seq(1,30)} 
    else if (k==2) {p <- seq(1,150)}
                      
    # Initialisation des critères d'évaluation
    sensibilite <- matrix(nrow = length(simulation)/2, ncol = length(methodes))
    fdr <- matrix(nrow = length(simulation)/2, ncol = length(methodes))
    stabilite <- matrix(nrow = ((length(simulation)/2)*(length(simulation)/2-1))/2, ncol = length(methodes))
    rmse <- matrix(nrow = length(simulation)/2, ncol = length(methodes))
    empower <- list(boruta_ep = matrix(0, length(simulation)/2, 5000), 
                    janitza_ep = matrix(0,  length(simulation)/2, 5000))
    time <- matrix(0, nrow = length(simulation)/2, ncol = length(methodes))
      
    vars_select <- list()
    RF <- list()
    
    i = 1
    
    # Application des méthodes sur la première moitié des réplicas
    for (replica in simulation[1:(length(simulation)/2)]){
      
      # Boruta 
      time[i,1] <- system.time({
        boruta <- Boruta(replica$x, replica$y)
        methodes$boruta_vs <- which(boruta$finalDecision=='Confirmed')
        modeles$boruta_rf <- randomForest(x = replica$x[,methodes$boruta_vs], y = replica$y)})[3]

      # Janitza
      time[i,2] <- system.time({
        PerVarImp1 <- CVPVI(replica$x, replica$y)
        janitza <- NTA(PerVarImp1$cv_varim)
        methodes$janitza_vs <- which(janitza$pvalue==0)
        modeles$janitza_rf <- randomForest(x = replica$x[,methodes$janitza_vs], y = replica$y)})[3]
      
      # # Altmann
      # time[i,3] <- system.time({
      #   rf <- randomForest(replica$x, replica$y, importance = T)
      #   PerVarImp2 <- PIMP(replica$x, replica$y, rForest = rf, S = 50)
      #   altmann <- PimpTest(PerVarImp2)
      #   methodes$altmann_vs <- which(altmann$pvalue==0)
      #   modeles$altmann_rf <- randomForest(x = replica$x[,methodes$altmann_vs], y = replica$y)})[3]
      # 
      # # VSURF
      # time[i,4] <- system.time({
      #   vsurf <- VSURF(replica$x, replica$y)
      #   methodes$vsurf_vs <- vsurf$varselect.interp
      #   modeles$vsurf_rf <- randomForest(x = replica$x[,methodes$vsurf_vs], y = replica$y)})[3]
      
      # CoV/VSURF
      # covsurf <- covsurf(replica$x, replica$y)
      # covsurf$vsurf_ptree$varselect.interp 
      
      # Sensibilité & FDR
      j <- 1
      for (vars in methodes){
        vp <- length(intersect(vars, p))
        fn <- length(p) - vp
        fp <- length(setdiff(vars, p))
        
        sensibilite[i,j] <- vp/(vp+fn)
        fdr[i,j] <- fp/(fp+vp)
        j <- j+1
      }
      
      vars_select[[i]] <- methodes
      RF[[i]] <- modeles
      i <- i+1
      pb$tick()
    }
    
    # Stabilité
    pair_index <- 1
    for (a in 1:(length(vars_select) - 1)) {
      for (b in (a + 1):length(vars_select)) {
        for (m in 1:length(methodes)) {
          vars_a <- vars_select[[a]][[m]]
          vars_b <- vars_select[[b]][[m]]
          intersection_ab <- length(intersect(vars_a, vars_b))
          union_ab <- length(union(vars_a, vars_b))
          stabilite[pair_index, m] <- intersection_ab / union_ab
        }
        pair_index <- pair_index + 1
      }
    }
    
    # RMSE 
    i = 1
    # Application des modèles séletcionnés sur l'autre moitié des réplicas
    for (replica in simulation[(length(simulation)/2+1):length(simulation)]){
      for (m in 1:length(methodes)){
        prediction <- predict(RF[[i]][[m]], newdata = replica$x[,vars_select[[i]][[m]]])
        rmse[i,m] <- sqrt(mean((replica$y-prediction)^2))
      }
      i <- i+1
    }
    
    # Empirical Power
    i <- 1
    for (v in 1:length(vars_select)){
      for (m in 1:length(methodes)){
        empower[[m]][i,vars_select[[i]][[m]]] <- empower[[m]][i,vars_select[[i]][[m]]] + 1
      }
      i <- i+1
    }
    
    # Moyennes des temps de calculs et de la puissance empirique pour un scénario
    for (e in 1:length(empower)){empower[[e]] <- colMeans(empower[[e]])}
    empower_tot[[k]] <- empower
    time_tot[k,] <- colMeans(time)
    
    # Enregistrement des valeurs pour calculer ensuite Q1, Q2, et Q3 pour les autres critères
    valeurs$sensibilite[[k]] <- sensibilite
    valeurs$fdr[[k]] <- fdr
    valeurs$stabilite[[k]] <- stabilite
    valeurs$rmse[[k]] <- rmse
    
    k <- k+1
    pb$tick()
  }
  
  # Calcul des médianes et intervalles interquartiles des critères restants
  criteres <- list(sensibilite = matrix(nrow = 2, ncol = length(methodes), dimnames = list(c(1,2), col_names)),
                     fdr = matrix(nrow = 2, ncol = length(methodes), dimnames = list(c(1,2), col_names)),
                     stabilite = matrix(nrow = 2, ncol = length(methodes), dimnames = list(c(1,2), col_names)),
                     rmse = matrix(nrow = 2, ncol = length(methodes), dimnames = list(c(1,2), col_names)))
  iq <- list(sensibilite = data.frame(), 
             fdr = data.frame(), 
             stabilite = data.frame(), 
             rmse = data.frame())
  for (c in 1:length(criteres)){
    for (i in 1:nrow(criteres[[c]])){
      for (j in 1:ncol(criteres[[c]])){
        quartiles <- round(quantile(valeurs[[c]][[i]][,j], probs = c(0.25,0.5,0.75)), 4)
        criteres[[c]][i,j] <- quartiles[2]
        iq[[c]][i,j] <- paste(quartiles[1], quartiles[3], sep = ';')
        colnames(iq[[c]])[j] <- col_names[j]
      }
    }
  }
  
  
  resultats <- list(sensibilité = list(valeurs = criteres[[1]], iq = iq$sensibilite), 
                    fdr = list(valeurs = criteres[[2]], iq = iq$fdr), 
                    stabilité = list(valeurs = criteres[[3]], iq = iq$stabilite), 
                    rmse = list(valeurs = criteres[[4]], iq = iq$rmse), 
                    empower = empower_tot,
                    time = time_tot)
  
  return(resultats)
}

graphiques <- function(resultats){
  
  Méthodes <- c('Boruta','Janitza')
  
  Q1 <- list(list(c(), c(), c(), c()),
             list(c(), c(), c(), c()))
  Q3 <- list(list(c(), c(), c(), c()),
             list(c(), c(), c(), c()))
  for (i in 1:length(resultats[1:4])){
    for (j in 1:nrow(resultats[[i]]$iq)){
      for (col in resultats[[i]]$iq[j,]){
        Q1[[j]][[i]] <- c(Q1[[j]][[i]], as.numeric(strsplit(col,';')[[1]][1]))
        Q3[[j]][[i]] <- c(Q3[[j]][[i]], as.numeric(strsplit(col,';')[[1]][2]))
      }
    }
  }
  
  # Sensibilité/FDR
  data11 <- data.frame(x = resultats$fdr$valeurs[1,],
                       y = resultats$sensibilité$valeurs[1,],
                       Q1_x = Q1[[1]][[2]],
                       Q3_x = Q3[[1]][[2]],
                       Q1_y = Q1[[1]][[1]],
                       Q3_y = Q3[[1]][[1]])
  
  graph11 <- ggplot(data11, aes(x = x, y = y, col = Méthodes)) + geom_point() +
    geom_errorbar(aes(xmin = Q1_x, xmax = Q3_x), width = 0.05, color = c('darkgreen','red')) +
    geom_errorbar(aes(ymin = Q1_y, ymax = Q3_y), width = 0.05, color = c('darkgreen','red')) +
    labs(title = "n = 10",
         x = "FDR",
         y = "Sensibilité") +
    coord_cartesian(xlim = c(0, 1), ylim = c(0,1)) +
    scale_color_manual(values = c("Boruta" = "darkgreen", "Janitza" = "red")) +
    scale_y_continuous(breaks = seq(0, 1, by = 0.2)) +
    theme_minimal()
  
  data12 <- data.frame(x = resultats$fdr$valeurs[2,],
                       y = resultats$sensibilité$valeurs[2,],
                       Q1_x = Q1[[2]][[2]],
                       Q3_x = Q3[[2]][[2]],
                       Q1_y = Q1[[2]][[1]],
                       Q3_y = Q3[[2]][[1]])
  
  graph12 <- ggplot(data12, aes(x = x, y = y, col = Méthodes)) + geom_point() +
    geom_errorbar(aes(xmin = Q1_x, xmax = Q3_x), width = 0.05, color = c('darkgreen','red')) +
    geom_errorbar(aes(ymin = Q1_y, ymax = Q3_y), width = 0.05, color = c('darkgreen','red')) +
    labs(title = "n = 50",
         x = "FDR",
         y = "Sensibilité") +
    coord_cartesian(xlim = c(0, 1), ylim = c(0,1)) +
    scale_color_manual(values = c("Boruta" = "darkgreen", "Janitza" = "red")) +
    scale_y_continuous(breaks = seq(0, 1, by = 0.2)) +
    theme_minimal()
  
  # Stabilité/RMSE
  data21 <- data.frame(x = resultats$rmse$valeurs[1,],
                       y = resultats$stabilité$valeurs[1,],
                       Q1_x = Q1[[1]][[4]],
                       Q3_x = Q3[[1]][[4]],
                       Q1_y = Q1[[1]][[3]],
                       Q3_y = Q3[[1]][[3]])
  graph21 <- ggplot(data21, aes(x = x, y = y, col = Méthodes)) + geom_point() +
    geom_errorbar(aes(xmin = Q1_x, xmax = Q3_x), width = 0.05, color = c('darkgreen','red')) +
    geom_errorbar(aes(ymin = Q1_y, ymax = Q3_y), width = 0.05, color = c('darkgreen', 'red')) +
    labs(title = "n = 10",
         x = "RMSE",
         y = "Stabilité") +
    coord_cartesian(xlim = c(1, 2), ylim = c(0,1)) +
    scale_color_manual(values = c("Boruta" = "darkgreen", "Janitza" = "red")) +
    scale_y_continuous(breaks = seq(0, 1, by = 0.2)) +
    theme_minimal()
  
  data22 <- data.frame(x = resultats$rmse$valeurs[2,],
                       y = resultats$stabilité$valeurs[2,],
                       Q1_x = Q1[[2]][[4]],
                       Q3_x = Q3[[2]][[4]],
                       Q1_y = Q1[[2]][[3]],
                       Q3_y = Q3[[2]][[3]])
  graph22 <- ggplot(data22, aes(x = x, y = y, col = Méthodes)) + geom_point() +
    geom_errorbar(aes(xmin = Q1_x, xmax = Q3_x), width = 0.05, color = c('darkgreen','red')) +
    geom_errorbar(aes(ymin = Q1_y, ymax = Q3_y), width = 0.05, color = c('darkgreen', 'red')) +
    labs(title = "n = 50",
         x = "RMSE",
         y = "Stabilité") +
    coord_cartesian(xlim = c(1, 2), ylim = c(0,1)) +
    scale_color_manual(values = c("Boruta" = "darkgreen", "Janitza" = "red")) +
    scale_y_continuous(breaks = seq(0, 1, by = 0.2)) +
    theme_minimal() 
  
  grid.arrange(graph11, graph12, graph21, graph22, ncol = 2)
  
  
  # Empirical Power
  data31 <- data.frame(x = seq(1,30),
                       Boruta = resultats$empower[[1]]$boruta_ep[seq(1,30)],
                       Janitza = resultats$empower[[1]]$janitza_ep[seq(1,30)])
  graph31 <- ggplot(data31, aes(x = x)) + 
    geom_line(aes(y = Boruta, color = "Boruta")) +
    geom_line(aes(y = Janitza, color = "Janitza")) +
    labs(title = 'n = 10', 
         x = "Variables d'intérêt",
         y = "Fréquence de sélection", 
         color = 'Méthodes') +
    scale_color_manual(values = c("Boruta" = "darkgreen", "Janitza" = "red"),
                       labels = c("Boruta", "Janitza")) +
    scale_y_continuous(breaks = seq(0, 1, by = 0.2)) +
    theme_minimal()
  
  data32 <- data.frame(x = seq(1,150),
                       Boruta = resultats$empower[[2]]$boruta_ep[seq(1,150)],
                       Janitza = resultats$empower[[2]]$janitza_ep[seq(1,150)])
  graph32 <- ggplot(data32, aes(x = x)) + 
    geom_line(aes(y = Boruta, color = "Boruta")) +
    geom_line(aes(y = Janitza, color = "Janitza")) +
    labs(title = 'n = 50',
         x = "Variables d'intérêt",
         y = "Fréquence de sélection", 
         color = 'Méthodes') +
    scale_color_manual(values = c("Boruta" = "darkgreen", "Janitza" = "red"),
                       labels = c("Boruta", "Janitza")) +
    scale_y_continuous(breaks = seq(0, 1, by = 0.2)) +
    theme_minimal()
  
  grid.arrange(graph31, graph32, ncol = 2)
}

# Test #
simulation <- simulation_regression(R=100)
resultats <- evaluation(simulation)
graphiques(resultats)

# Détails des résultats #
resultats$sensibilité
resultats$fdr
resultats$stabilité
resultats$rmse
resultats$time/60



# # Test des méthodes individuellement 

# ## BORUTA ##
# test0 <- Boruta(simu[[2]][[3]]$x, simu[[2]][[3]]$y)
# test1 <- which(test0$finalDecision=='Confirmed')
# length(intersect(test1, p2))
# 
# ## JANITZA ##
# test1 <- CVPVI(simu[[2]][[2]]$x, simu[[2]][[2]]$y) 
# test2 <- NTA(test1$cv_varim)
# # Visualisation
# mean(test2$pvalue)
# median(test2$pvalue)
# plot(sort(test2$pvalue)[1:5000])
# # Valeur seuil
# which(test2$pvalue==0)
# which(test2$pvalue<median(test2$pvalue))
# 
## ALTMANN ##
# RF <- randomForest(simu[[1]][[5]]$x, simu[[1]][[5]]$y, importance = TRUE)
# test3 <- PIMP(simu[[1]][[5]]$x, simu[[1]][[5]]$y, RF)
# test4 <- PimpTest(test3)
# which(test4$pvalue==0)
# 
# # VSURF
# test5 <- VSURF(simu[[2]][[1]]$x, simu[[2]][[1]]$y)
# test5$varselect.thres
# sort(test5$varselect.interp)
# test5$mean.perf
# 
# COVSURF
# test6 <- covsurf(simu[[1]][[5]]$x, simu[[1]][[5]]$y, kval = c(2:10))
# sort(test6$vsurf_ptree$varselect.interp)
