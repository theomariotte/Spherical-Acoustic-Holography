%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% validation du calcul des coefficients de Fourier par résolution d'un
% système linéaire mal conditionné
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear; clc;
close all;

%% Paramètres
data_path = 'data/';
% Tracé de la pression acoustique
% plt_typ : ('real') real part ; ('dB') decibels ; ('abs') absolute value
plt_typ = 'dB';

% Méthode de simulation
sim_method = 'Propagator';

% Microphone array radius
a = 15e-2;
% fréquence de travail [Hz]
f = 721;
% Maximum order (Bessel, Hankel, SH)
Nmax = 4;
% Noise amp
nAmp = 0.00;
% color limits
clims = [70 100];

%%% Autres paramètres
% débit de la source [m^3/s]
Q = 1e-4;
% masse volumique du milieu
rho_air = 1.2;
% célérité du milieu 
c = 340;

%% Source location

% cartesian coordinates [m]
xs = [.3];
ys = [.3];
zs = [.3];

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
                'ReconstructRadius',a,...
                'c',c,...
                'rho',rho_air,...
                'Q',Q,...
                'method',sim_method,...
                'incidentOnly',0,...
                'doplot',0);
[P] = generateSimu(Rm,Rs,pp_simu,0); 

% RSB computation
num_RSB = max(max(abs(P)));
RSB = 20*log10(num_RSB/nAmp);

fprintf('RSB = %2.0f dB\n',RSB)

%% Calcul des coefficients de Fourier
idx = 1;
Ymat = zeros(Nmic,(Nmax+1)^2);
Pmn_th = zeros((Nmax+1)^2,1);
k = 2*pi*f / c;

for n = 0 : Nmax
    
    % pour Fourier analytique
    [hn_r0,~,~] = SphericalHankel1(n,k*rhos);
    [~,dhn_a,~] = SphericalHankel1(n,k*a);
    
   for m = -n : n
       % Matrice des harmoniques sphériques
        Y_tmp = getSphericalHarmonics(theta,phi,n,m); 
        Ymat(:,idx) = Y_tmp;
                
        % coefficients de Fourier analytiques
        Y_source = getSphericalHarmonics(thetas,phis,n,m);
%         Pmn_th(idx) = -4*pi*hn_r0/((k*a)^2*dhn_a) * conj(Y_source); 
        Pmn_th(idx) = -4*pi/(k*a)^2 * ( hn_r0 / dhn_a )* conj(Y_source); 
        % update
        idx = idx + 1;        
   end
end

P_meas = P + nAmp * randn(size(P));
Pmn = sphericalTFSolver(Ymat,P_meas);

%% Figures

% pression simulée
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

% cofficients de Fourier
hh = figure('Name','Fourier coefficient');
hold on
plot(20*log10(abs(Pmn_th)),...
    'kx','linewidth',2,...
    'markersize',7);
plot(20*log10(abs(Pmn)),...
    'o','linewidth',2,...
    'markersize',7,...
    'color',[.6 .6 .6]);    
hold off
xlabel('Coefficient index')
ylabel('20log_{10}(|P_{mn}|)')
title(sprintf('Fourier coefficients : ka = %2.1f',k*a));
legend('Theoretical','Computed')
grid on
set(gca,'ylim',[-100 20],'fontsize',12)

fig_path = 'Figures/';
fname = sprintf('FOURIER_f_%d_RSB_%d',f,round(RSB));
printFigFmt(hh,fig_path,fname,'eps');

