spherical acoustic holography/plot

Fonctions permettant de tracer diverses grandeurs. 

Contenu du dossier :

dessin_13HP.m (JH Thomas) : tracé de la source de bruit

getOutFileName.m : obtention d'un nom de fichier (vérifie s'il existe, si oui, un nouveau nom est automatiquement généré ; utile pour les figures). 

pressureMeasurementVisu : visualisation du champ de pression mesuré par l'antenne sur une sphère. Le champ est interpolé sur un nombre de points déterminé par l'utilisateur puis tracÈ en couleurs
sur une sphère.

printFigFmt.m : enregistre une figure à un format donné. 

RadialFuncVisu.m : visualisation des fonctions radiales (Bessel Hankel sphériques). Cette fonction
permet de tracer les fonctions et leurs dérivées à partir d'un vecteur r. 

SHvisualization : tracÈ des harmoniques sphériques. 3 types de visualisation sont possibles : tracÈ
en 3D, projection sur la surface d'une sphère, projection selon les trois plans de l'espace (2D).
La fonction prend en argument une grille d'anges et une matrice contenant les valeurs d'un harmonique sphérique associées à la grille.