
source('Fonctions.R')


# 1 - Simulation
simulations <- simulation('regression', R = 10)

# 2 - Évaluation
resultats <- evaluation(simulations)
# Détails
resultats$sensibilite
resultats$fdr
resultats$stabilite
resultats$erreur.pred
resultats$times.methodes
resultats$time

# 3 - Graphiques
graphiques(resultats)