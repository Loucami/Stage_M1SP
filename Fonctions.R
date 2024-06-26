

# 1 - Fonction de simulation des données 
simulation <- function(type, N = 100, P = 5000, p = 6, M = c(10,50), R = 100){
  
  if (!require('mlbench', character.only = TRUE, quietly = TRUE)) {
    install.packages('mlbench')
    require('mlbench', character.only = TRUE, quietly = TRUE)
  }
  
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


# 2 - Fonction d'évaluation
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
    col_names <- c('Boruta', 'Vita', 'Altmann', 'VSURF')
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
                       Vita_vs = c(),
                       altmann_vs = c(),
                       vsurf_vs = c())
      
      modeles <- list(boruta_rf = c(),
                      Vita_rf = c(),
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
                      Vita_ep = matrix(0,  length(simulation)/2, 5000),
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
        
        # Vita
        time[i,2] <- system.time({
          PerVarImp1 <- CVPVI(replica$x, replica$y, parallel = TRUE)
          Vita <- NTA(PerVarImp1$cv_varim)
          methodes$Vita_vs <- which(Vita$pvalue==0)
          if (length(methodes$Vita_vs)!=0){modeles$Vita_rf <- randomForest(x = replica$x[, methodes$Vita_vs, drop=FALSE], y = replica$y)}})[3]
        
        # Altmann
        time[i,3] <- system.time({
          rf <- randomForest(replica$x, replica$y, importance = T)
          PerVarImp2 <- PIMP(replica$x, replica$y, rForest = rf, S = 50, parallel = TRUE)
          altmann <- PimpTest(PerVarImp2, para = F)
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
  time.function <- paste(minutes, "min", secondes, "s")
  
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


# 3 - Fonction d'évalution avec parallélisation de Boruta
evaluation_parallel <- function(simulations){
  
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
      total = length(simulations)*(length(simulations[[1]])/2)+2*length(simulations)+1
    )
    
    # Critères d'évaluation
    col_names <- c('Boruta', 'Vita', 'Altmann', 'VSURF')
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
                       Vita_vs = c(),
                       altmann_vs = c(),
                       vsurf_vs = c())
      
      modeles <- list(boruta_rf = c(),
                      Vita_rf = c(),
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
                      Vita_ep = matrix(0,  length(simulation)/2, 5000),
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
        
        # Vita
        time[i,2] <- system.time({
          PerVarImp1 <- CVPVI(replica$x, replica$y, parallel = TRUE)
          Vita <- NTA(PerVarImp1$cv_varim)
          methodes$Vita_vs <- which(Vita$pvalue==0)
          if (length(methodes$Vita_vs)!=0){modeles$Vita_rf <- randomForest(x = replica$x[, methodes$Vita_vs, drop=FALSE], y = replica$y)}})[3]
        
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
  })[3]
  
  minutes <- floor(time.function/60)
  secondes <- round((time.function/60-minutes)*60)
  time.function <- paste(minutes, "min", secondes, "s")
  
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


# 4 - Fonction permettant la création de graphiques
graphiques <- function(resultats){
  
  packages <- c("ggplot2", 'gridExtra')
  
  for (pkg in packages) {
    if (!require(pkg, character.only = TRUE, quietly = TRUE)) {
      install.packages(pkg)
      require(pkg, character.only = TRUE, quietly = TRUE)
    }
  }
  
  Méthodes <- c('Boruta', 'Vita', 'Altmann', 'VSURF')
  
  Q1 <- list(list(c(), c(), c(), c(), c()),
             list(c(), c(), c(), c(), c()))
  Q3 <- list(list(c(), c(), c(), c(), c()),
             list(c(), c(), c(), c(), c()))
  
  for (i in 1:length(resultats[1:5])){
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
    geom_errorbar(aes(xmin = Q1_x, xmax = Q3_x), width = 0.05, color = c('darkgreen', 'red', 'darkblue', 'orange')) +
    geom_errorbar(aes(ymin = Q1_y, ymax = Q3_y), width = 0.05, color = c('darkgreen', 'red', 'darkblue', 'orange')) +
    labs(title = "Groupes n = 10",
         x = "TFD",
         y = "Sensibilité") +
    coord_cartesian(xlim = c(0, 1), ylim = c(0,1)) +
    scale_color_manual(values = c('Altmann' = 'darkblue', "Boruta" = "darkgreen", "Vita" = "red", 'VSURF' = 'orange')) +
    scale_y_continuous(breaks = seq(0, 1, by = 0.2)) +
    theme_minimal() +
    theme(plot.title = element_text(hjust = 0.5))
  
  data12 <- data.frame(x = resultats$fdr$valeurs[2,],
                       y = resultats$sensibilité$valeurs[2,],
                       Q1_x = Q1[[2]][[2]],
                       Q3_x = Q3[[2]][[2]],
                       Q1_y = Q1[[2]][[1]],
                       Q3_y = Q3[[2]][[1]])
  
  graph12 <- ggplot(data12, aes(x = x, y = y, col = Méthodes)) + geom_point() +
    geom_errorbar(aes(xmin = Q1_x, xmax = Q3_x), width = 0.05, color = c('darkgreen', 'red', 'darkblue', 'orange')) +
    geom_errorbar(aes(ymin = Q1_y, ymax = Q3_y), width = 0.05, color = c('darkgreen', 'red', 'darkblue', 'orange')) +
    labs(title = "Groupes n = 50",
         x = "TFD",
         y = "Sensibilité") +
    coord_cartesian(xlim = c(0, 1), ylim = c(0,1)) +
    scale_color_manual(values = c('Altmann' = 'darkblue', "Boruta" = "darkgreen", "Vita" = "red", 'VSURF' = 'orange')) +
    scale_y_continuous(breaks = seq(0, 1, by = 0.2)) +
    theme_minimal() +
    theme(plot.title = element_text(hjust = 0.5))
  
  # Stabilité/Erreur de prédiction
  if (resultats$type == 'regression'){xbornes <- c(1,2)}
  else if (resultats$type == 'classification'){xbornes <- c(0,1)}
  
  data21 <- data.frame(x = resultats$erreur.pred$valeurs[1,],
                       y = resultats$stabilité$valeurs[1,],
                       Q1_x = Q1[[1]][[4]],
                       Q3_x = Q3[[1]][[4]],
                       Q1_y = Q1[[1]][[3]],
                       Q3_y = Q3[[1]][[3]])
  graph21 <- ggplot(data21, aes(x = x, y = y, col = Méthodes)) + geom_point() +
    geom_errorbar(aes(xmin = Q1_x, xmax = Q3_x), width = 0.05, color = c('darkgreen', 'red', 'darkblue', 'orange')) +
    geom_errorbar(aes(ymin = Q1_y, ymax = Q3_y), width = 0.05, color = c('darkgreen', 'red', 'darkblue', 'orange')) +
    labs(x = "EC",
         y = "Stabilité") +
    coord_cartesian(xlim = xbornes, ylim = c(0,1)) +
    scale_color_manual(values = c('Altmann' = 'darkblue', "Boruta" = "darkgreen", "Vita" = "red", 'VSURF' = 'orange')) +
    scale_y_continuous(breaks = seq(0, 1, by = 0.2)) +
    theme_minimal() 
  
  data22 <- data.frame(x = resultats$erreur.pred$valeurs[2,],
                       y = resultats$stabilité$valeurs[2,],
                       Q1_x = Q1[[2]][[4]],
                       Q3_x = Q3[[2]][[4]],
                       Q1_y = Q1[[2]][[3]],
                       Q3_y = Q3[[2]][[3]])
  graph22 <- ggplot(data22, aes(x = x, y = y, col = Méthodes)) + geom_point() +
    geom_errorbar(aes(xmin = Q1_x, xmax = Q3_x), width = 0.05, color = c('darkgreen', 'red', 'darkblue', 'orange')) +
    geom_errorbar(aes(ymin = Q1_y, ymax = Q3_y), width = 0.05, color = c('darkgreen', 'red', 'darkblue', 'orange')) +
    labs(x = "EC",
         y = "Stabilité") +
    coord_cartesian(xlim = xbornes, ylim = c(0,1)) +
    scale_color_manual(values = c('Altmann' = 'darkblue', "Boruta" = "darkgreen", "Vita" = "red", 'VSURF' = 'orange')) +
    scale_y_continuous(breaks = seq(0, 1, by = 0.2)) +
    theme_minimal()
  
  grid.arrange(graph11, graph12, graph21, graph22, ncol = 2)
  
  if (nrow(resultats$nb.fp$valeurs)>2){
    data.null1 <- data.frame(x = resultats$erreur.pred$valeurs[3,],
                             y = resultats$nb.fp$valeurs[3,],
                             Q1_x = Q1[[3]][[4]],
                             Q3_x = Q3[[3]][[4]],
                             Q1_y = Q1[[3]][[5]],
                             Q3_y = Q3[[3]][[5]])
    graph.null1 <- ggplot(data.null1, aes(x = x, y = y, col = Méthodes)) + geom_point() +
      geom_errorbar(aes(xmin = Q1_x, xmax = Q3_x), width = 0.05, color = c('darkgreen', 'red', 'darkblue', 'orange')) +
      geom_errorbar(aes(ymin = Q1_y, ymax = Q3_y), width = 0.05, color = c('darkgreen', 'red', 'darkblue', 'orange')) +
      labs(title = "n = 10",
           x = "Erreur de prédiction",
           y = "Nombre de faux positifs") +
      coord_cartesian(xlim = c(0,2)) +
      scale_color_manual(values = c('Altmann' = 'darkblue', "Boruta" = "darkgreen", "Vita" = "red", 'VSURF' = 'orange')) +
      theme_minimal() +
      theme(plot.title = element_text(hjust = 0.5))
    
    data.null2 <- data.frame(x = resultats$erreur.pred$valeurs[4,],
                             y = resultats$nb.fp$valeurs[4,],
                             Q1_x = Q1[[4]][[4]],
                             Q3_x = Q3[[4]][[4]],
                             Q1_y = Q1[[4]][[5]],
                             Q3_y = Q3[[4]][[5]])
    graph.null2 <- ggplot(data.null2, aes(x = x, y = y, col = Méthodes)) + geom_point() +
      geom_errorbar(aes(xmin = Q1_x, xmax = Q3_x), width = 0.05, color = c('darkgreen', 'red', 'darkblue', 'orange')) +
      geom_errorbar(aes(ymin = Q1_y, ymax = Q3_y), width = 0.05, color = c('darkgreen', 'red', 'darkblue', 'orange')) +
      labs(title = "n = 50",
           x = "Erreur de prédiction",
           y = "Nombre de faux positifs") +
      coord_cartesian(xlim = c(0,2)) +
      scale_color_manual(values = c('Altmann' = 'darkblue', "Boruta" = "darkgreen", "Vita" = "red", 'VSURF' = 'orange')) +
      theme_minimal() +
      theme(plot.title = element_text(hjust = 0.5))
    
    grid.arrange(graph.null1, graph.null2, ncol = 2)
  }
  
  
  # Empirical Power
  data31 <- data.frame(x = seq(1,30),
                       Altmann = resultats$empower[[1]]$altmann_ep[seq(1,30)],
                       Boruta = resultats$empower[[1]]$boruta_ep[seq(1,30)],
                       Vita = resultats$empower[[1]]$janitza_ep[seq(1,30)],
                       VSURF = resultats$empower[[1]]$vsurf_ep[seq(1,30)])
  graph31 <- ggplot(data31, aes(x = x)) + 
    geom_line(aes(y = Altmann, color = "Altmann")) +
    geom_line(aes(y = Boruta, color = "Boruta")) +
    geom_line(aes(y = Vita, color = "Vita")) +
    geom_line(aes(y = VSURF, color = "VSURF")) +
    labs(title = 'Groupes n = 10', 
         x = "Variables d'intérêt",
         y = "Fréquence de sélection", 
         color = 'Méthodes') +
    scale_color_manual(values = c('Altmann' = 'darkblue', "Boruta" = "darkgreen", "Vita" = "red", 'VSURF' = 'orange'),
                       labels = c('Altmann', 'Boruta', "Vita", 'VSURF')) +
    scale_y_continuous(breaks = seq(0, 1, by = 0.2)) +
    theme_minimal() +
    theme(plot.title = element_text(hjust = 0.5))
  
  data32 <- data.frame(x = seq(1,150),
                       Altmann = resultats$empower[[2]]$altmann_ep[seq(1,150)],
                       Boruta = resultats$empower[[2]]$boruta_ep[seq(1,150)],
                       Vita = resultats$empower[[2]]$janitza_ep[seq(1,150)],
                       VSURF = resultats$empower[[2]]$vsurf_ep[seq(1,150)])
  graph32 <- ggplot(data32, aes(x = x)) + 
    geom_line(aes(y = Altmann, color = "Altmann")) +
    geom_line(aes(y = Boruta, color = "Boruta")) +
    geom_line(aes(y = Vita, color = "Vita")) +
    geom_line(aes(y = VSURF, color = "VSURF")) +
    labs(title = 'Groupes n = 50',
         x = "Variables d'intérêt",
         y = "Fréquence de sélection", 
         color = 'Méthodes') +
    scale_color_manual(values = c('Altmann' = 'darkblue', "Boruta" = "darkgreen", "Vita" = "red", 'VSURF' = 'orange'),
                       labels = c('Altmann', 'Boruta', "Vita", 'VSURF')) +
    scale_y_continuous(breaks = seq(0, 1, by = 0.2)) +
    theme_minimal() +
    theme(plot.title = element_text(hjust = 0.5))
  
  grid.arrange(graph31, graph32, ncol = 2)
}

