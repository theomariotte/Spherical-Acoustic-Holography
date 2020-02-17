Spherical acoustic holography/_tests

Dossier contenant les scripts de test de diff�rentes partie du code de reconstruction. 

test_validation_Fourier.m : calcul des coefficients de fourier (pour la m�thode d'holographie sph�rique uniquement, pas pour les sources �quivalentes)

test_jn_hn.m : test de calcul des fonctions de Bessel et de Hankel sph�riques + d�riv�es

SESM_from_meas.m : calcul de la reconstruction � partir de mesures (impl�mentation en cours)

test_simu_meas.m : calcul d'une simulation ( = calculer le champ sur la surface d'une sph�re rigide 
en des points discrets lorsqu'un ou plusieurs monopoles rayonnent). Perte de g�n�rer son propre 
r�seau de microphones ou d'utiliser les coordonn�es de la sph�re du LAUM. 

test_spherical_harmonics : calcul des harmoniques sph�riques et validation avec Williams

test_SESM : calcul d'une reconstruction avec un ou plusieurs monop�le(s) avec plusieurs m�thodes de r�gularisation (tikhonov ou lsqr)

test_optimalRegularization : calcul du param�tre de reg. Optimal en utilisant Generalized Cross Validation


