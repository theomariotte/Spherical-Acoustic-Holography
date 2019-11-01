Spherical acoustic holography/_tests

Dossier contenant les scripts de test de différentes partie du code de reconstruction. 

test_Fourier_coeff.m : test de calcul des coefficients de Fourier sphériques (TF sphérique). 
Le calcul réalisé dans cette fonction n'a pas été validé car les résultats de la simulation
ne sont pas encore vérifiés. 

test_jn_hn.m : test de calcul des fonctions de Bessel et de Hankel sphériques + dérivées

test_propagator.m : calcul du propagateur donné par Williams dans "Intensity vector reconstruction ..."
Deux calculs sont réalisés et permettent de valider le calcul des fonctions de Bessel et Hankel 
sphériques ainsi que le calcul des harmoniques sphériques. 

test_simu_meas.m : calcul d'une simulation ( = calculer le champ sur la surface d'une sphère rigide 
en des points discrets lorsqu'un ou plusieurs monopoles rayonnent). Perte de générer son propre 
réseau de microphones ou d'utiliser les coordonnées de la sphère du LAUM. (Résultats étranges,
à vérifier)

test_spherical_harmonics : est du calcul des harmoniques sphériques, comparaison avec quelques 
valeurs analytiques données par Williams et tracé. 


