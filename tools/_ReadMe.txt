spherical acoustic holography / tools

Fonctions utiles pour la m�thode des sources �quivalentes sph�riques et autre.

Contenu du dossier :

cartesianCoordinates.m : passage en coordonn�es cart�siennes depuis le sph�rique

conversion_principale.m : code de conversion des donn�es exp�rimentales

convert_data_TMDS.m : utile � la conversion des donn�es exp�rimentales �galement

equivalentSourcesGrid.m : construction d'une grille de sources �quivalentes (non test�)

fieldInterpolation : interpolation d'un champ sur une nouvelle grille (principalement pour 
faciliter le trac� du champ)

getFreeFieldFRF.m : fonction de transfert en champ libre (fonction de Green) pour la reconstruction du champ de pression. 

getGreenNeumannFRF.m : fonction de transfert (fonction de Green Neumann) pour r�soudre le probl�me inverse (calcul du d�bit des sources �quivalentes)

getSphericalHarmonics.m : fonction calculant un harmonique sph�rique pour un degr� et un ordre donn�. Elle prend en argument une grille d'anges (azimut et �l�vation). 

randomESlocation.m : r�partition des sources �quivalentes al�atoirement sur une sph�re de rayon donn� (fonctionne mais ne peut pas �tre utilis� dans l'algorithme de reconstruction pour l'instant car matrice 3D)

solveIllPosedProblem.m : r�solution d'un syst�me mal pos�/conditionn�. R�gularisation de Tikhonov + LSQR

SphericalBessel.m : fonction calculant les fonctions de Bessel sph�riques de types 1 et 2.

SphericalHanel1.m : calcul des fonctions de Hankel du premier type. 

SphericalHanel2.m : calcul des fonctions de Hankel du second type. 

trace_sphere_36.m : trac� de la position des microphones

visu_data_TMDS.m : visualisation des donn�es exp�rimentales