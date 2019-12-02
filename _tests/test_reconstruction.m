% test simulation mesures
clear; clc;
close all;

%% Paramètres

data_path = 'data/';
% Tracé de la pression acoustique
% plt_typ : ('real') real part ; ('dB') decibels ; ('abs') absolute value
plt_typ = 'abs';

% Méthode de simulation
% ('GN') : utilisation de la fonction de Green Neumann (calcul possible
% uniquement sur la frontière) ; 
% ('Propagator') Utilise le propagateur
% défini par Williams dans 'Intensity vector reconstruction'. 
sim_method = 'Propagator';

% Microphone array radius
a = 15e-2;
% fréquence de travail [Hz]
f = 100;
% Maximum order (Bessel, Hankel, SH)
Nmax = 6;
% Noise amp
nAmp = 0;

%%% Autres paramètres
% débit de la source [m^3/s]
Q = 1e-3;
% masse volumique du milieu
rho_air = 1.2;
% célérité du milieu 
c = 340;

%% Source location

% cartesian coordinates [m]
xs = 0.5;
ys = 0.5;
zs = 0.5;

Rs = [xs ys zs];
[rhos,phis,thetas] = sphericalCoordinates(xs,ys,zs);

%% Load microphone locations

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
        Pmn_th(idx) = -4*pi*hn_r0/((k*a)^2*dhn_a) * conj(Y_source); 
        
        % update
        idx = idx + 1;        
   end
end

pp_solve = struct('regularization',1,...
           'alpha',1e-10,...
           'doplot',0,...
           'compute_condition',1);
P_meas = P + nAmp * randn(size(P));
[Pmn,cond] = solveIllPosedProblem(Ymat,P,pp_solve);


%% Affichage 

% Coefficients de Fourier
hh = figure('Name','Fourier coefficient');
hold on
    plot(20*log10(abs(Pmn_th)),'kx');
    plot(20*log10(abs(Pmn)),'ro');    
hold off
xlabel('n')
ylabel('20log_{10}(|P_{mn}|)')
title(sprintf('Fourier coefficients : ka = %2.1f',k*a));
legend('Theoretical','Computed')
grid on
set(gca,'ylim',[-100 20])







pp_interp = struct('sphereRadius',a,...
            'numAngle',30,...
            'interpType','natural',...
            'doplot',0);
        
pp_plot = struct('showSource',1,...
                'showMic',0,...
                'plt_typ',plt_typ,...
                'clim',[],...
                'fontSize',12,...
                'freq',f);     
            
handle = pressureMeasurementVisu(P_meas,Rm,Rs,pp_interp,pp_plot);


% Pressure field on a sphere with data interpolation (does not work yet)
% cartesian frame


% [Xsph,Ysph,Zsph] = sph2cart(phiq,thetaq,P2plot);
% 
% % tst = convhul(Xsph,Ysph);
% 
% % draw a unit sphere
% [xs,ys,zs]=sphere(max(size(P2plot))-1);
% xs = a * xs;
% ys = a * ys;
% zs = a * zs;
% 
% h2 = figure('Name','Pressure on the sphere');
% hold on
% surf(xs,ys,zs,P2plot);
% colormap('jet')
% shading('interp')
% % plot3(Rm(:,1),Rm(:,2),Rm(:,3),'ro','color',[.9 .5 .2],'linewidth',3);
% plot3(Rs(:,1),Rs(:,2),Rs(:,3),'ro','linewidth',3)
% hold off
% colorbar
% xlabel('X','interpreter','latex')
% ylabel('Y','interpreter','latex')
% zlabel('Z','interpreter','latex')
% axis('equal')

% h = figure('Name','Set up');
% hold on
% [x,y,z] = sphere(300);
% x = (a-0.005) * x;
% y = (a-0.005) * y;
% z = (a-0.005) * z;
% 
% surf(x,y,z,ones(size(x)));
% colormap('gray')
% shading('interp')
% camlight
% lighting phong
% plot3(Rm(:,1),Rm(:,2),Rm(:,3),'x','color',[.9 .5 .2],'linewidth',3);
% plot3(Rs(:,1),Rs(:,2),Rs(:,3),'ro','linewidth',3)
% grid on
% axis equal
% set(gca,'visible','off');
% hold off
% xlabel('X','interpreter','latex')
% ylabel('Y','interpreter','latex')
% zlabel('Z','interpreter','latex')
% axis('equal')
