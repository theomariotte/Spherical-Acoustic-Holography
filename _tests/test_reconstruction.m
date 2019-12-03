% test simulation mesures
clear; clc;
close all;

%% Paramètres

% rayon de la sphère de reconstruction
r_reconstruct = 30e-2;

data_path = 'data/';
% Tracé de la pression acoustique
% plt_typ : ('real') real part ; ('dB') decibels ; ('abs') absolute value
plt_typ = 'dB';

% Méthode de simulation
% ('GN') : utilisation de la fonction de Green Neumann (calcul possible
% uniquement sur la frontière) ; 
% ('Propagator') Utilise le propagateur
% défini par Williams dans 'Intensity vector reconstruction'. 
sim_method = 'GN';

% Microphone array radius
a = 15e-2;
% fréquence de travail [Hz]
f = 100;
% Maximum order (Bessel, Hankel, SH)
Nmax = 3;
% Noise amp
nAmp = 0.01;
% color limits
clims = [70 100];

%%% Autres paramètres
% débit de la source [m^3/s]
Q = 1e-3;
% masse volumique du milieu
rho_air = 1.2;
% célérité du milieu 
c = 340;

%% Source location

% cartesian coordinates [m]
xs = [0.5];
ys = [0.2];
zs = [0.5];

Rs = [xs ys zs];
[rhos,phis,thetas] = sphericalCoordinates(xs,ys,zs);

%% Chargement ds coordonnées des microphones (cartésien)

Nmic = 36;
filename =  '3Dcam36Officiel.txt';
fID = fopen([data_path filename],'r');
sizeA = [Nmic 3];
mic_loc_tmp = fscanf(fID,'%f');
fclose(fID);

Rm = transpose(reshape(mic_loc_tmp,[3 Nmic]));  

[rho,phi,theta] = sphericalCoordinates(Rm(:,1),Rm(:,2),Rm(:,3));

%% calculate the soundfield at each microphone

pp_simu = struct('MaxOrder',Nmax,...
                'freq',f,...
                'SphereRadius',a,...
                'c',c,...
                'rho',rho_air,...
                'Q',Q,...
                'method',sim_method,...
                'doplot',0);
          
[P] = generateSimu(Rm,Rs,pp_simu);        

pp_simu = struct('MaxOrder',Nmax,...
                'freq',f,...
                'SphereRadius',r_reconstruct,...
                'c',c,...
                'rho',rho_air,...
                'Q',Q,...
                'method',sim_method,...
                'doplot',0);

[P_cmp] = generateSimu(Rm,Rs,pp_simu);        

P_meas = P + nAmp * randn(size(P));
%% reconstruction du champ acoustique 

pp_reconstruction = struct('a',a,...
                        'R',r_reconstruct,...
                        'freq',f,...
                        'c',c,...
                        'maxOrder',Nmax,...
                        'incidentOnly',0);
    
[reconstruct_SF] = SNAH(theta,phi,P_meas,pp_reconstruction);

%% Figures

% pression simulée sur la sphère (interpolée pour une lecture plus simple)
pp_interp = struct('sphereRadius',a,...
            'numAngle',30,...
            'interpType','natural',...
            'doplot',0);
        
pp_plot = struct('showSource',1,...
                'showMic',0,...
                'plt_typ',plt_typ,...
                'clim',clims,...
                'fontSize',12,...
                'freq',f);     
            
handle1 = pressureMeasurementVisu(P_meas,Rm,Rs,pp_interp,pp_plot);

% Pression reconstruite par méthode d'holographie acoustique sphérique
pp_interp = struct('sphereRadius',r_reconstruct,...
            'numAngle',30,...
            'interpType','natural',...
            'doplot',0);
        
pp_plot = struct('showSource',1,...
                'showMic',0,...
                'plt_typ',plt_typ,...
                'clim',clims,...
                'fontSize',12,...
                'freq',f);     
            
handle2 = pressureMeasurementVisu(reconstruct_SF,Rm,Rs,pp_interp,pp_plot);
title('Recontructed Sound Field');

% Pression reconstruite attendue (calculée par simulation)
pp_interp = struct('sphereRadius',r_reconstruct,...
            'numAngle',30,...
            'interpType','natural',...
            'doplot',0);
        
pp_plot = struct('showSource',1,...
                'showMic',0,...
                'plt_typ',plt_typ,...
                'clim',clims,...
                'fontSize',12,...
                'freq',f);     
            
handle3 = pressureMeasurementVisu(P_cmp,Rm,Rs,pp_interp,pp_plot);
title('Expected Recontructed Sound Field');

