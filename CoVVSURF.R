
source('Fonctions.R')
library(mlbench)
library(CoVVSURF)
library(ClustOfVar)
library(mlbench)

resultats <- list(regression = list(),
                  classification = list())

for (type in c('regression', 'classification')){
  
  res.sensibilite <- matrix(nrow = 2, ncol = 3, dimnames = list(c('S1','S2'), c('Q1','Médiane','Q3')))
  res.fdr <- matrix(nrow = 2, ncol = 3, dimnames = list(c('S1','S2'), c('Q1','Médiane','Q3')))
  res.replicas <- list(scenario1 = list(), scenario2 = list())
  
  simulations <- simulation(type, R = 1)
  
  for (i in 1:2){
    
    if (i == 1) {vars <- paste0('X',seq(1,30))} else {vars <- paste0('X',seq(1,150))}
    
    simu <- simulations$data[[i]][[1]]
    sensibilite <- c()
    fdr <- c()
    
    for (j in 1:2){
      
      clusters <- kmeansvar(simu$x, init = 100)
      tree <- hclustvar(clusters$scores)
      selection <- covsurf(clusters$scores, simu$y, tree = tree, ncores = 6)
      
      cluselect <- list()
      varselect <- c()
      
      k <- 1
      
      for (cluster in selection$vsel){
        
        print(cluster)
        cluselect[[k]] <- clusters$var[[paste0('cluster', cluster)]]
        print(varselect)
        varselect <- c(varselect, names(which(clusters$var[[paste0('cluster', cluster)]][,1]>0.5)))
        k <- k+1
      } 
      
      sensibilite <- c(sensibilite, length(intersect(varselect, vars)) / length(vars))
      fdr <- c(fdr, length(setdiff(varselect, vars)) / length(varselect))
      
      res.replicas[[i]][[j]] <- list(kopt = selection$kopt,
                                     varselect = varselect,
                                     cluselect = cluselect)
    }
    
    res.sensibilite[i,] <- round(quantile(sensibilite, probs = c(0.25,0.5,0.75), na.rm = TRUE), 3)
    res.fdr[i,] <- round(quantile(fdr, probs = c(0.25,0.5,0.75), na.rm = TRUE), 3)
  }
  
  resultats[[type]] <- list(sensibilité = res.sensibilite,
                            fdr = res.fdr,
                            détails = res.replicas)
}


save(resultats, file = 'resultatsHCLUSTVAR.RData')

