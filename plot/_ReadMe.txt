spherical acoustic holography/plot

Fonctions permettant de tracer diverses grandeurs. 

Contenu du dossier :

pressureMeasurementVisu : visualisation du champ de pression mesuré par l'antenne sur une sphère. 
Le champ est interpolé sur un nombre de points déterminé par l'utilisateur puis tracé en couleurs
sur une sphère.

RadialFuncVisu.m : visualisation des fonctions radiales (Bessel Hankel sphériques). Cette fonction
permet de tracer les fonctions et leurs dérivées à partir d'un vecteur r. 

SHvisualization : tracé des harmoniques sphériques. 3 types de visulaisation sont possibles : tracé
en 3D, prjection sur la surface d'une sphère, projection selon les trois plans de l'espace (2D).
La fonction prend en argument une grille d'anges et une matrice contenant les valeurs d'un harmonique
sphérique associés à la grille.