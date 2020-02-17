% test calcul coefficients de Fourier

clear; clc; 
close all

tmp_home = '/Users/theomariotte/Documents/01_work/ENSIM/5A/Projet/Spherical-Acoustic-Holography/';
data_path = [tmp_home 'data/'];

%% paramètres
sim_method = 'GN';
a = 15e-2;
% fréquence de travail [Hz]
f = 1000;
% Maximum order (Bessel, Hankel, SH)
Nmax = 4;

%%% Autres paramètres
% débit de la source [m^3/s]
Q = 1e-3;
% masse volumique du milieu
rho = 1;
% célérité du milieu 
c = 340;

% position source
xs = 0;
ys = 0;
zs = 0.5;




%% chargement des positions des micros
Nmic = 36;
filename =  '3Dcam36Officiel.txt';
fID = fopen([data_path filename],'r');
sizeA = [Nmic 3];
mic_loc_tmp = fscanf(fID,'%f');
fclose(fID);

Rm = transpose(reshape(mic_loc_tmp,[3 Nmic]));  
dointerp = 1;



%% Calcul de la pression simulée

pp_simu = struct('MaxOrder',Nmax,...
               'freq',f,...
               'SphereRadius',a,...
               'ReconstructRadius',a,...
               'c',c,...
               'rho',rho,...
               'Q',Q,...
               'method',sim_method,...
               'incidentOnly',0,...
               'doplot',0);
Rs = [xs ys zs];            
[P] = generateSimu(Rm,Rs,pp_simu);    

%% calcul des harmoniques sphériques

[az,elev,r] = cart2sph(Rm(:,1),Rm(:,2),Rm(:,3));

idx_SH = 1;

for n = 0 : Nmax    
    for m = -n : n
        for imic = 1 : Nmic
            Ymn(imic,idx_SH) = getSphericalHarmonics(elev(imic),az(imic),n,m);
        end
        idx_SH = idx_SH+1;
    end    
end

%% Calcul des coefficients de Fourier

[Pmn,cond] = solveIllPosedProblem(Ymn,P);

