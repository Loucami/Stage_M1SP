
# Simulations
library(mlbench)
# Méthodes
library(randomForest)
library(VSURF)
# Fonctions
source('Fonctions.R')


simulations <- simulation('classification', R = 50)


# Mtry
# i <- 1
# mtry <- c(70,seq(250,5000,250))
# resultats <- matrix(nrow = 2, ncol = length(mtry), dimnames = list(c(1,2), c(70,seq(250,5000,250))))
# for (simu in simulations$data){
#   j <- 1
#   for (m in mtry){
#     OOBerror <- c()
#     for (replica in simu){
#       res <- randomForest(replica$x, replica$y, mtry = m)
#       OOBerror <- c(OOBerror, mean(res$err.rate[,1]))
#       pb$tick()
#     }
#     resultats[i,j] <- mean(OOBerror)
#     j <- j+1
#   }
#   i <- i+1
# }
load('Résultats/resultatsMtry.RData')
plot(y = resultats[1,], x = c(70,100,150,200,seq(250,5000,250)), xlab = 'Mtry', ylab = 'Erreur OOB', type = 'b')
plot(y = resultats[2,], x = c(70,100,150,200,seq(250,5000,250)), xlab = 'Mtry', ylab = 'Erreur OOB', type = 'b')


# NSD 
# i <- 1
# paramètres <- c(1,0.8,0.6,0.4,0.2,0)
# resultats <- list(scenario1 = matrix(nrow = length(paramètres), ncol = 3, dimnames = list(c(), c('Q1','Médiane','Q3'))),
#                   scenario2 = matrix(nrow = length(paramètres), ncol = 3, dimnames = list(c(), c('Q1','Médiane','Q3'))))
# for (simu in simulations$data){
#   if (i==1){vars <- 30} else {vars <- 150}
#   j <- 1
#   for (para in paramètres){
#     sensibilité <- c()
#     for (replica in simu) {
#       vsurf <- VSURF(replica$x, replica$y, parallel = TRUE, nsd = para)
#       varselect <- vsurf$varselect.interp
#       sensibilité <- c(sensibilité, length(intersect(varselect, seq(1,vars)))/vars)
#     }
#     resultats[[i]][j,] <- round(quantile(sensibilité, probs = c(0.25,0.5,0.75), na.rm = TRUE), 3)
#     j <- j+1
#   }
#   i <- i+1
# }
load('Résultats/resultatsSD.RData')
plot(y = resultats$scenario2[,2], x = seq(1,0,-0.1), xlab = 'SD', ylab = 'Sensibilité', type = 'b')


# Comparaison des méthodes avec Mtry = 1000 et NSD = 0 pour VSURF
load('Résultats/resultatsC1000-0.RData')
graphiques(resultats)

