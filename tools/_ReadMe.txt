spherical acoustic holography / tools

Différents outils utile à l'application de la méthode de l'holographie acoustique sphérique. 

Contenu du dossier :

cartesianCoordinates.m : passage en coordonnées cartésiennes depuis le sphérique

fieldInterpolation : interpolation d'un champ sur une nouvelle grille (principalement pour 
faciliter le tracé du champ)

SNAH.m : Spherical Acoustic Holography solver (not implemented yet!)

getsSphericalHarmonics.m : fonction calculant un hamonique sphérique pour un degré et un ordre
donné. Elle prend en argument une grille d'anges (azimut et élévation). 

solveIllPosedProblem.m : fonction permettant de résoudre un problème mal posé (pseudo inversion de 
matrice et régularisation de Tychonov). < non validée ! >

SphericalBessel.m : fonction calculant les fonctions de Bessel sphériques de types 1 et 2.

SphericalHanel1.m : calcul des fonctions de Hankel du premier type. 