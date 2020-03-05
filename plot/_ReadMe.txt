spherical acoustic holography/plot

Fonctions permettant de tracer diverses grandeurs. 

Contenu du dossier :

dessin_13HP.m (JH Thomas) : trac� de la source de bruit

getOutFileName.m : obtention d'un nom de fichier (v�rifie s'il existe, si oui, un nouveau nom est automatiquement g�n�r� ; utile pour les figures). 

pressureMeasurementVisu : visualisation du champ de pression mesur� par l'antenne sur une sph�re. Le champ est interpol� sur un nombre de points d�termin� par l'utilisateur puis trac� en couleurs
sur une sph�re.

printFigFmt.m : enregistre une figure � un format donn�. 

RadialFuncVisu.m : visualisation des fonctions radiales (Bessel Hankel sph�riques). Cette fonction
permet de tracer les fonctions et leurs d�riv�es � partir d'un vecteur r. 

SHvisualization : trac� des harmoniques sph�riques. 3 types de visualisation sont possibles : trac�
en 3D, projection sur la surface d'une sph�re, projection selon les trois plans de l'espace (2D).
La fonction prend en argument une grille d'anges et une matrice contenant les valeurs d'un harmonique sph�rique associ�es � la grille.