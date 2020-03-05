spherical acoustic holography / tools

Fonctions utiles pour la méthode des sources équivalentes sphériques et autre.

Contenu du dossier :

cartesianCoordinates.m : passage en coordonnées cartésiennes depuis le sphérique

conversion_principale.m : code de conversion des données expérimentales

convert_data_TMDS.m : utile à la conversion des données expérimentales également

equivalentSourcesGrid.m : construction d'une grille de sources équivalentes (non testé)

fieldInterpolation : interpolation d'un champ sur une nouvelle grille (principalement pour 
faciliter le tracé du champ)

getFreeFieldFRF.m : fonction de transfert en champ libre (fonction de Green) pour la reconstruction du champ de pression. 

getGreenNeumannFRF.m : fonction de transfert (fonction de Green Neumann) pour résoudre le problème inverse (calcul du débit des sources équivalentes)

getSphericalHarmonics.m : fonction calculant un harmonique sphérique pour un degré et un ordre donné. Elle prend en argument une grille d'anges (azimut et élévation). 

randomESlocation.m : répartition des sources équivalentes aléatoirement sur une sphère de rayon donné (fonctionne mais ne peut pas être utilisé dans l'algorithme de reconstruction pour l'instant car matrice 3D)

solveIllPosedProblem.m : résolution d'un système mal posé/conditionné. Régularisation de Tikhonov + LSQR

SphericalBessel.m : fonction calculant les fonctions de Bessel sphériques de types 1 et 2.

SphericalHanel1.m : calcul des fonctions de Hankel du premier type. 

SphericalHanel2.m : calcul des fonctions de Hankel du second type. 

trace_sphere_36.m : tracé de la position des microphones

visu_data_TMDS.m : visualisation des données expérimentales