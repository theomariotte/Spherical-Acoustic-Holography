spherical acoustic holography / tools

Diff�rents outils utile � l'application de la m�thode de l'holographie acoustique sph�rique. 

Contenu du dossier :

cartesianCoordinates.m : passage en coordonn�es cart�siennes depuis le sph�rique

fieldInterpolation : interpolation d'un champ sur une nouvelle grille (principalement pour 
faciliter le trac� du champ)

SNAH.m : Spherical Acoustic Holography solver (not implemented yet!)

getsSphericalHarmonics.m : fonction calculant un hamonique sph�rique pour un degr� et un ordre
donn�. Elle prend en argument une grille d'anges (azimut et �l�vation). 

solveIllPosedProblem.m : fonction permettant de r�soudre un probl�me mal pos� (pseudo inversion de 
matrice et r�gularisation de Tychonov). < non valid�e ! >

SphericalBessel.m : fonction calculant les fonctions de Bessel sph�riques de types 1 et 2.

SphericalHanel1.m : calcul des fonctions de Hankel du premier type. 