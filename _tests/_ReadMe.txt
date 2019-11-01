Spherical acoustic holography/_tests

Dossier contenant les scripts de test de diff�rentes partie du code de reconstruction. 

test_Fourier_coeff.m : test de calcul des coefficients de Fourier sph�riques (TF sph�rique). 
Le calcul r�alis� dans cette fonction n'a pas �t� valid� car les r�sultats de la simulation
ne sont pas encore v�rifi�s. 

test_jn_hn.m : test de calcul des fonctions de Bessel et de Hankel sph�riques + d�riv�es

test_propagator.m : calcul du propagateur donn� par Williams dans "Intensity vector reconstruction ..."
Deux calculs sont r�alis�s et permettent de valider le calcul des fonctions de Bessel et Hankel 
sph�riques ainsi que le calcul des harmoniques sph�riques. 

test_simu_meas.m : calcul d'une simulation ( = calculer le champ sur la surface d'une sph�re rigide 
en des points discrets lorsqu'un ou plusieurs monopoles rayonnent). Perte de g�n�rer son propre 
r�seau de microphones ou d'utiliser les coordonn�es de la sph�re du LAUM. (R�sultats �tranges,
� v�rifier)

test_spherical_harmonics : est du calcul des harmoniques sph�riques, comparaison avec quelques 
valeurs analytiques donn�es par Williams et trac�. 


