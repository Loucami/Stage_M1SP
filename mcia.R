
install.packages('randomForest', repos="https://cloud.r-project.org/")
install.packages('VSURF', repos="https://cloud.r-project.org/")
install.packages('Boruta', repos="https://cloud.r-project.org/")
install.packages('vita', repos="https://cloud.r-project.org/")
install.packages('mlbench', repos="https://cloud.r-project.org/")
install.packages('ggplot2', repos="https://cloud.r-project.org/")
install.packages('gridExtra', repos="https://cloud.r-project.org/")
install.packages('progress', repos="https://cloud.r-project.org/")
library(randomForest)
library(VSURF)
library(Boruta)
library(vita)
library(mlbench)
library(ggplot2)
library(gridExtra)
library(progress)
library(foreach)
library(doParallel)


simulation <- function(type, N = 100, P = 5000, p = 6, M = c(10,50), R = 100){
  
  simulations = list()
  k = 1
  
  # Scénarios dépendants
  for (m in M){
    data = list()
    
    for (r in 1:R){
      # On génère x et Y 
      if (type == 'regression'){
        x <- matrix(runif(p*N), ncol = p)
        epsilon <- rnorm(N, 0, 0.2) 
        Y <- 0.25*exp(4*x[,1]) + 4/(1+exp(-20*(x[,2]-0.5))) + 3*x[,3] + epsilon
      } else if (type == 'classification') {
        mlb.data <- mlbench.threenorm(N, d = 3)
        x <- cbind(mlb.data$x, matrix(runif(N*(p-3),-3,3), nrow = N))
        Y <- mlb.data$classes
      }
      
      # On crée p*m variables corrélées
      vars <- matrix(NA, N, p*m)
      if (type == 'regression') {sd <- 0.3} else {sd <- sqrt(2)}
      for (i in 1:p) {
        for (j in 1:m){
          delta <- rnorm(N, 0, sd) 
          vars[,(i-1)*m+j] <- x[,i] + (0.01 + (0.5 * (j - 1)) / (m - 1))*delta
        } 
      }
      
      # On crée notre échantillon X en agrégeant les variables corrélées (vars) avec des variables indépendantes (vars_sup)
      if (type == 'regression'){vars_sup <- matrix(runif(N*(P-p*m)), ncol = P-p*m)} 
      else {vars_sup <- matrix(runif(N*(P-p*m),-3,3), ncol = P-p*m)}
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
  #     if (type == 'regression'){
  #       x_ind <- matrix(rnorm(3*N, 0, 0.2), ncol = 3)
  #       epsilon <- rnorm(N, 0, 0.2)
  #       Y <- 0.25*exp(4*x_ind[,1]) + 4/(1+exp(-20*(x_ind[,2]-0.5))) + 3*x_ind[,3] + epsilon
  #     } else {
  #       mlb.data <- mlbench.threenorm(N, d = 3)
  #       Y <- mlb.data$classes
  #     }
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
  
  return(list(data = simulations, type = type))
}
evaluation <- function(simulations){
  
  time.function <- system.time({
    
    packages <- c('randomForest', "Boruta", "vita", 'VSURF', 'progress')
    
    for (pkg in packages) {
      if (!require(pkg, character.only = TRUE, quietly = TRUE)) {
        install.packages(pkg)
        require(pkg, character.only = TRUE, quietly = TRUE)
      }
    }
    
    type <- simulations$type
    simulations <- simulations$data
    
    pb <- progress_bar$new(
      format = "[:bar] :percent ETA: :eta",
      total = length(simulations)*(length(simulations[[1]])/2)+length(simulations)+1
    )
    
    # Critères d'évaluation
    col_names <- c('Boruta', 'Janitza', 'Altmann', 'VSURF')
    if(length(simulations)==2){row_names <- c(1,2)} else {row_names <- c(1,2,3,4)}
    empower_tot <- list()
    time_tot <- matrix(nrow = length(simulations), ncol = length(col_names), dimnames = list(row_names, col_names))
    valeurs <- list(sensibilite = list(), 
                    fdr = list(),
                    stabilite = list(),
                    erreur.pred = list(),
                    nb.fp = list(),
                    f1score = list())
    
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
      if (k==1) {p <- seq(1,30)} 
      else if (k==2) {p <- seq(1,150)}
      
      # Initialisation des critères d'évaluation
      sensibilite <- matrix(nrow = length(simulation)/2, ncol = length(methodes))
      fdr <- matrix(nrow = length(simulation)/2, ncol = length(methodes))
      stabilite <- matrix(nrow = ((length(simulation)/2)*(length(simulation)/2-1))/2, ncol = length(methodes))
      erreur.pred <- matrix(nrow = length(simulation)/2, ncol = length(methodes))
      empower <- list(boruta_ep = matrix(0, length(simulation)/2, 5000), 
                      janitza_ep = matrix(0,  length(simulation)/2, 5000),
                      altmann_ep = matrix(0,  length(simulation)/2, 5000),
                      vsurf_ep = matrix(0,  length(simulation)/2, 5000))
      time <- matrix(0, nrow = length(simulation)/2, ncol = length(methodes))
      
      nb.fp <- matrix(nrow = length(simulation)/2, ncol = length(methodes))
      f1score <- matrix(nrow = length(simulation)/2, ncol = length(methodes))
      
      vars_select <- list()
      RF <- list()
      
      i = 1
      
      # Application des méthodes sur la première moitié des réplicas
      for (replica in simulation[1:(length(simulation)/2)]){
        
        # Boruta 
        time[i,1] <- system.time({
          boruta <- Boruta(replica$x, replica$y)
          boruta <- TentativeRoughFix(boruta)
          methodes$boruta_vs <- which(boruta$finalDecision=='Confirmed')
          modeles$boruta_rf <- randomForest(x = replica$x[,methodes$boruta_vs], y = replica$y)})[3]
        
        # Janitza
        time[i,2] <- system.time({
          PerVarImp1 <- CVPVI(replica$x, replica$y, parallel = TRUE)
          janitza <- NTA(PerVarImp1$cv_varim)
          methodes$janitza_vs <- which(janitza$pvalue==0)
          if (length(methodes$janitza_vs)!=0){modeles$janitza_rf <- randomForest(x = replica$x[, methodes$janitza_vs, drop=FALSE], y = replica$y)}})[3]
        
        # Altmann
        time[i,3] <- system.time({
          rf <- randomForest(replica$x, replica$y, importance = T)
          PerVarImp2 <- PIMP(replica$x, replica$y, rForest = rf, S = 50, parallel = TRUE)
          altmann <- PimpTest(PerVarImp2, para = T)
          methodes$altmann_vs <- which(altmann$pvalue==0)
          if (length(methodes$altmann_vs)!=0){modeles$altmann_rf <- randomForest(x = replica$x[, methodes$altmann_vs, drop=FALSE], y = replica$y)}})[3]
        
        # VSURF
        time[i,4] <- system.time({
          vsurf <- VSURF(replica$x, replica$y, parallel = TRUE)
          methodes$vsurf_vs <- vsurf$varselect.interp
          modeles$vsurf_rf <- randomForest(x = replica$x[,methodes$vsurf_vs], y = replica$y)})[3]
        
        # Sensibilité & FDR
        j <- 1
        for (vars in methodes){
          vp <- length(intersect(vars, p))
          fn <- length(p) - vp
          fp <- length(setdiff(vars, p))
          
          nb.fp[i,j] <- fp
          f1score[i,j] <- vp/(vp+(1/2)*(fn+fp))
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
      
      # Erreur de prédiction 
      i = 1
      # Application des modèles séletcionnés sur l'autre moitié des réplicas
      for (replica in simulation[(length(simulation)/2+1):length(simulation)]){
        for (m in 1:length(methodes)){
          if (length(vars_select[[i]][[m]])!=0){
            prediction <- predict(RF[[i]][[m]], newdata = replica$x[, vars_select[[i]][[m]], drop=FALSE])
            if (type == 'regression'){erreur.pred[i,m] <- sqrt(mean((replica$y-prediction)^2))} # RMSE
            else if (type == 'classification'){erreur.pred[i,m] <- sum(replica$y!=prediction)/length(prediction)} # Erreur de classification
          }
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
      valeurs$erreur.pred[[k]] <- erreur.pred
      valeurs$nb.fp[[k]] <- nb.fp
      valeurs$f1score[[k]] <- f1score
      
      k <- k+1
      pb$tick()
    }
    
    # Calcul des médianes et intervalles interquartiles des critères restants
    criteres <- list(sensibilite = matrix(nrow = length(simulations), ncol = length(methodes), dimnames = list(row_names, col_names)),
                     fdr = matrix(nrow = length(simulations), ncol = length(methodes), dimnames = list(row_names, col_names)),
                     stabilite = matrix(nrow = length(simulations), ncol = length(methodes), dimnames = list(row_names, col_names)),
                     erreur.pred = matrix(nrow = length(simulations), ncol = length(methodes), dimnames = list(row_names, col_names)),
                     nb.fp = matrix(nrow = length(simulations), ncol = length(methodes), dimnames = list(row_names, col_names)),
                     f1score = matrix(nrow = length(simulations), ncol = length(methodes), dimnames = list(row_names, col_names)))
    
    iq <- list(sensibilite = data.frame(),
               fdr = data.frame(), 
               stabilite = data.frame(), 
               erreur.pred = data.frame(),
               nb.fp = data.frame(),
               f1score = data.frame())
    
    for (c in 1:length(criteres)){
      for (i in 1:nrow(criteres[[c]])){
        for (j in 1:ncol(criteres[[c]])){
          quartiles <- round(quantile(valeurs[[c]][[i]][,j], probs = c(0.25,0.5,0.75), na.rm = TRUE), 3)
          criteres[[c]][i,j] <- quartiles[2]
          iq[[c]][i,j] <- paste(quartiles[1], quartiles[3], sep = ';')
          colnames(iq[[c]])[j] <- col_names[j]
        }
      }
    }
  })[3]
  
  minutes <- floor(time.function/60)
  secondes <- round((time.function/60-minutes)*60)
  time.function <- paste(minutes, "min", secondes, "s\n")
  
  resultats <- list(sensibilité = list(valeurs = criteres[[1]], iq = iq$sensibilite), 
                    fdr = list(valeurs = criteres[[2]], iq = iq$fdr), 
                    stabilité = list(valeurs = criteres[[3]], iq = iq$stabilite), 
                    erreur.pred = list(valeurs = criteres[[4]], iq = iq$erreur.pred), 
                    nb.fp = list(valeurs = criteres[[5]], iq = iq$nb.fp),
                    f1score = list(valeurs = criteres[[6]], iq = iq$f1score),
                    empower = empower_tot,
                    times.methodes = time_tot,
                    time = time.function,
                    type = type)
  
  return(resultats)
}
evaluation_parallel <- function(simulations){
  
  type <- simulations$type
  simulations <- simulations$data
  
  pb <- progress_bar$new(
    format = "[:bar] :percent ETA: :eta",
    total = length(simulations)*(length(simulations[[1]])/2)+2*length(simulations)+1
  )
  
  # Critères d'évaluation
  col_names <- c('Boruta', 'Janitza', 'Altmann', 'VSURF')
  if(length(simulations)==2){row_names <- c(1,2)} else {row_names <- c(1,2,3,4)}
  empower_tot <- list()
  time_tot <- matrix(nrow = length(simulations), ncol = length(col_names), dimnames = list(row_names, col_names))
  valeurs <- list(sensibilite = list(), 
                  fdr = list(),
                  stabilite = list(),
                  erreur.pred = list(),
                  nb.fp = list(),
                  f1score = list())
  
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
    if (k==1) {p <- seq(1,30)} 
    else if (k==2) {p <- seq(1,150)}
    
    # Initialisation des critères d'évaluation
    sensibilite <- matrix(nrow = length(simulation)/2, ncol = length(methodes))
    fdr <- matrix(nrow = length(simulation)/2, ncol = length(methodes))
    stabilite <- matrix(nrow = ((length(simulation)/2)*(length(simulation)/2-1))/2, ncol = length(methodes))
    erreur.pred <- matrix(nrow = length(simulation)/2, ncol = length(methodes))
    empower <- list(boruta_ep = matrix(0, length(simulation)/2, 5000), 
                    janitza_ep = matrix(0,  length(simulation)/2, 5000),
                    altmann_ep = matrix(0,  length(simulation)/2, 5000),
                    vsurf_ep = matrix(0,  length(simulation)/2, 5000))
    time <- matrix(0, nrow = length(simulation)/2, ncol = length(methodes))
    
    nb.fp <- matrix(nrow = length(simulation)/2, ncol = length(methodes))
    f1score <- matrix(nrow = length(simulation)/2, ncol = length(methodes))
    
    vars_select <- list()
    RF <- list()
    
    # Initialisation du cluster pour la parallélisation 
    num_cores <- detectCores() - 1
    cl <- makeCluster(num_cores)
    registerDoParallel(cl)
    
    # Parallélisation de Boruta sur la première moitié des réplicas
    valeurs_boruta <- foreach(replica = simulation[1:(length(simulation)/2)], .packages = c("Boruta", "randomForest")) %dopar% {
      
      # Boruta 
      time.B <- system.time({
        boruta <- Boruta(replica$x, replica$y)
        boruta <- TentativeRoughFix(boruta)
        vars <- which(boruta$finalDecision=='Confirmed')
        rf <- randomForest(x = replica$x[,vars], y = replica$y)})[3]
      
      # Sensibilité & FDR 
      vp <- length(intersect(vars, p))
      fn <- length(p) - vp
      fp <- length(setdiff(vars, p))
      
      nb.fp.B <- fp
      f1score.B <- vp/(vp+(1/2)*(fn+fp))
      sensibilite.B <- vp/(vp+fn)
      fdr.B <- fp/(fp+vp)
      
      list(vars = vars, rf = rf, time = time.B, nb.fp = nb.fp.B, f1score = f1score.B, sensibilite = sensibilite.B, fdr = fdr.B)
    }
    
    stopCluster(cl)
    pb$tick()
    
    i = 1
    
    # Application des autres méthodes sur la première moitié des réplicas
    for (replica in simulation[1:(length(simulation)/2)]){
      
      # Janitza
      time[i,2] <- system.time({
        PerVarImp1 <- CVPVI(replica$x, replica$y, parallel = TRUE)
        janitza <- NTA(PerVarImp1$cv_varim)
        methodes$janitza_vs <- which(janitza$pvalue==0)
        if (length(methodes$janitza_vs)!=0){modeles$janitza_rf <- randomForest(x = replica$x[, methodes$janitza_vs, drop=FALSE], y = replica$y)}})[3]
      
      # Altmann
      time[i,3] <- system.time({
        rf <- randomForest(replica$x, replica$y, importance = T)
        PerVarImp2 <- PIMP(replica$x, replica$y, rForest = rf, S = 50, parallel = TRUE)
        altmann <- PimpTest(PerVarImp2, para = T)
        methodes$altmann_vs <- which(altmann$pvalue==0)
        if (length(methodes$altmann_vs)!=0){modeles$altmann_rf <- randomForest(x = replica$x[, methodes$altmann_vs, drop=FALSE], y = replica$y)}})[3]
      
      # VSURF
      time[i,4] <- system.time({
        vsurf <- VSURF(replica$x, replica$y, parallel = TRUE)
        methodes$vsurf_vs <- vsurf$varselect.interp
        modeles$vsurf_rf <- randomForest(x = replica$x[,methodes$vsurf_vs], y = replica$y)})[3]
      
      # Sensibilité & FDR
      j <- 2
      for (vars in methodes[2:length(methodes)]){
        vp <- length(intersect(vars, p))
        fn <- length(p) - vp
        fp <- length(setdiff(vars, p))
        
        nb.fp[i,j] <- fp
        f1score[i,j] <- vp/(vp+(1/2)*(fn+fp))
        sensibilite[i,j] <- vp/(vp+fn)
        fdr[i,j] <- fp/(fp+vp)
        j <- j+1
      }
      
      vars_select[[i]] <- methodes
      RF[[i]] <- modeles
      i <- i+1
      pb$tick()
    }
    
    # Enregistrement des données de Boruta
    for (r in seq_along(simulation[1:(length(simulation)/2)])){
      vars_select[[r]]$boruta_vs <- valeurs_boruta[[r]]$vars
      RF[[r]]$boruta_rf <- valeurs_boruta[[r]]$rf 
      time[r,1] <- valeurs_boruta[[r]]$time
      nb.fp[r,1] <- valeurs_boruta[[r]]$nb.fp
      f1score[r,1] <- valeurs_boruta[[r]]$f1score
      sensibilite[r,1] <- valeurs_boruta[[r]]$sensibilite
      fdr[r,1] <- valeurs_boruta[[r]]$fdr
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
    
    # Erreur de prédiction 
    i = 1
    # Application des modèles séletcionnés sur l'autre moitié des réplicas
    for (replica in simulation[(length(simulation)/2+1):length(simulation)]){
      for (m in 1:length(methodes)){
        if (length(vars_select[[i]][[m]])!=0){
          prediction <- predict(RF[[i]][[m]], newdata = replica$x[, vars_select[[i]][[m]], drop=FALSE])
          if (type == 'regression'){erreur.pred[i,m] <- sqrt(mean((replica$y-prediction)^2))} # RMSE
          else if (type == 'classification'){erreur.pred[i,m] <- sum(replica$y!=prediction)/length(prediction)} # Erreur de classification
        }
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
    valeurs$erreur.pred[[k]] <- erreur.pred
    valeurs$nb.fp[[k]] <- nb.fp
    valeurs$f1score[[k]] <- f1score
    
    k <- k+1
    pb$tick()
  }
  
  # Calcul des médianes et intervalles interquartiles des critères restants
  criteres <- list(sensibilite = matrix(nrow = length(simulations), ncol = length(methodes), dimnames = list(row_names, col_names)),
                   fdr = matrix(nrow = length(simulations), ncol = length(methodes), dimnames = list(row_names, col_names)),
                   stabilite = matrix(nrow = length(simulations), ncol = length(methodes), dimnames = list(row_names, col_names)),
                   erreur.pred = matrix(nrow = length(simulations), ncol = length(methodes), dimnames = list(row_names, col_names)),
                   nb.fp = matrix(nrow = length(simulations), ncol = length(methodes), dimnames = list(row_names, col_names)),
                   f1score = matrix(nrow = length(simulations), ncol = length(methodes), dimnames = list(row_names, col_names)))
  
  iq <- list(sensibilite = data.frame(),
             fdr = data.frame(), 
             stabilite = data.frame(), 
             erreur.pred = data.frame(),
             nb.fp = data.frame(),
             f1score = data.frame())
  
  for (c in 1:length(criteres)){
    for (i in 1:nrow(criteres[[c]])){
      for (j in 1:ncol(criteres[[c]])){
        quartiles <- round(quantile(valeurs[[c]][[i]][,j], probs = c(0.25,0.5,0.75), na.rm = TRUE), 3)
        criteres[[c]][i,j] <- quartiles[2]
        iq[[c]][i,j] <- paste(quartiles[1], quartiles[3], sep = ';')
        colnames(iq[[c]])[j] <- col_names[j]
      }
    }
  }
  
  resultats <- list(sensibilité = list(valeurs = criteres[[1]], iq = iq$sensibilite), 
                    fdr = list(valeurs = criteres[[2]], iq = iq$fdr), 
                    stabilité = list(valeurs = criteres[[3]], iq = iq$stabilite), 
                    erreur.pred = list(valeurs = criteres[[4]], iq = iq$erreur.pred), 
                    nb.fp = list(valeurs = criteres[[5]], iq = iq$nb.fp),
                    f1score = list(valeurs = criteres[[6]], iq = iq$f1score),
                    empower = empower_tot,
                    time = time_tot,
                    type = type)
  
  return(resultats)
}


simulations <- simulation('regression', R = 4)

resultats <- evaluation(simulations)
resultats$times.methodes
resultats$time
resultats$sensibilité

# temps.parallel <- system.time(resultats <- evaluation_parallel(simulations))
# resultats$time
# resultats$sensibilité
# resultats$stabilité
# minutes <- floor(temps.parallel[3]/60)
# secondes <- (temps.parallel[3]/60-minutes)*60
# cat(minutes, "min", secondes, "secondes\n")

simulations <- simulations$data
temps <- system.time({
  RF <- randomForest(simulations[[2]][[1]]$x, simulations[[2]][[1]]$y, importance = TRUE)
  resultats <- PIMP(simulations[[2]][[1]]$x, simulations[[2]][[1]]$y, RF, S = 200, parallel = TRUE)
  resultats <- PimpTest(resultats, para = F)
})
temps
length(which(resultats$pvalue==0)) 
which(resultats$pvalue==0)