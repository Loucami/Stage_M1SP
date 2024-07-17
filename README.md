
## Stage M1SP : Comparaison de méthodes de sélection de variables basées sur des forêts aléatoires.

Ce dépôt Github renseigne les différents codes utilisés lors de mon stage de Master 1 en Santé Publique, réalisé d'avril à juillet 2024. Il vise à détailler au maximum les étapes de comparaison des méthodes de sélection de variables suivantes : Boruta, Altmann, Vita et VSURF. 


### Fonctions.R

Fichier contenant les fonctions créées afin de réaliser la comparaison. 


### Résultats.R

Fichier d'application des méthodes de sélection. 


### CoVVSURF.R

Fichier tentant une implémentation de la méthode COV/VSURF en plusieurs étapes, afin d'obtenir des temps de calculs convenables. Dans un premier temps, un clustering est réalisé sur les variables afin de réduire considérablement les dimensions des données, puis un arbre hiérarchique ascendant est construit. Finalement, cette arbre est utilisé par la fonction CoVVSURF pour réaliser la sélection. 


### Résultats 

Dossier contenant les différents résultats obtenus.


### Graphiques

Dossier contenant les graphiques illustrant les résultats obtenues.


