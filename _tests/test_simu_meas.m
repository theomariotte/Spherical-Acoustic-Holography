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
sim_method = 'GN';

% Microphone array radius
a = 15e-2;
% fréquence de travail [Hz]
f = 500;
% Maximum order (Bessel, Hankel, SH)
Nmax = 6;

%%% Autres paramètres
% débit de la source [m^3/s]
Q = 1e-3;
% masse volumique du milieu
rho = 1.2;
% célérité du milieu 
c = 340;

%% Source location

% cartesian coordinates [m]
xs = [ .3 ; 0 ; -.6];
ys = [ .3 ; .5 ; 0];
zs = [ .3 ; -.3 ; -.3];

Rs = [xs ys zs];

%% Load microphone locations
Nmic = 36;
filename =  '3Dcam36Officiel.txt';
fID = fopen([data_path filename],'r');
sizeA = [Nmic 3];
mic_loc_tmp = fscanf(fID,'%f');
fclose(fID);

Rm = transpose(reshape(mic_loc_tmp,[3 Nmic]));  

%% calculate the soundfield at each microphone

pp_simu = struct('MaxOrder',Nmax,...
                'freq',f,...
                'SphereRadius',a,...
                'c',c,...
                'rho',rho,...
                'Q',Q,...
                'method',sim_method,...
                'doplot',0);

            
[P] = generateSimu(Rm,Rs,pp_simu);        

%% Affichage 
pp_interp = struct('sphereRadius',a,...
            'numAngle',30,...
            'interpType','natural',...
            'doplot',1);
        
pp_plot = struct('showSource',0,...
                'showMic',0,...
                'plt_typ',plt_typ,...
                'clim',[],...
                'fontSize',12,...
                'freq',f);     
            
handle = pressureMeasurementVisu(P,Rm,Rs,pp_interp,pp_plot);


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
