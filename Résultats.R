
source('Fonctions.R')
load(file = 'Résultats/regression.RData')
load(file = 'Résultats/classification.RData')
load(file = 'Résultats/resultatsC1000-0.RData')

# 1 - Simulation
simulations <- simulation('regression', R = 10)

# 2 - Évaluation
resultats <- evaluation(simulations)
# Détails
resultats$sensibilité
resultats$fdr
resultats$stabilité
resultats$erreur.pred
resultats$empower
resultats$f1score
resultats$times.methodes
resultats$time

# 3 - Graphiques
graphiques(resultats)
