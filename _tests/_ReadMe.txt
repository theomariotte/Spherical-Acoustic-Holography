Spherical acoustic holography/_tests

Dossier contenant les scripts de test de différentes partie du code de reconstruction. 

test_validation_Fourier.m : calcul des coefficients de fourier (pour la méthode d'holographie sphérique uniquement, pas pour les sources équivalentes)

test_jn_hn.m : test de calcul des fonctions de Bessel et de Hankel sphériques + dérivées

SESM_from_meas.m : calcul de la reconstruction à partir de mesures (implémentation en cours)

test_simu_meas.m : calcul d'une simulation ( = calculer le champ sur la surface d'une sphËre rigide 
en des points discrets lorsqu'un ou plusieurs monopoles rayonnent). Perte de gÈnÈrer son propre 
rÈseau de microphones ou d'utiliser les coordonnÈes de la sphËre du LAUM. 

test_spherical_harmonics : calcul des harmoniques sphériques et validation avec Williams

test_SESM : calcul d'une reconstruction avec un ou plusieurs monopôle(s) avec plusieurs méthodes de régularisation (tikhonov ou lsqr)

test_optimalRegularization : calcul du paramètre de reg. Optimal en utilisant Generalized Cross Validation


