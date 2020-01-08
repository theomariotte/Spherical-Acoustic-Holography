% test simulation mesures
clear; clc;
close all;

%% Param�tres

data_path = 'data/';
% Trac� de la pression acoustique
% plt_typ : ('real') real part ; ('dB') decibels ; ('abs') absolute value
plt_typ = 'dB';

% M�thode de simulation
% ('GN') : utilisation de la fonction de Green Neumann (calcul possible
% uniquement sur la fronti�re) ; 
% ('Propagator') Utilise le propagateur
% d�fini par Williams dans 'Intensity vector reconstruction'. 
sim_method = 'GN';

% Microphone array radius
a = 15e-2;
% fr�quence de travail [Hz]
f = 721;
% Maximum order (Bessel, Hankel, SH)
Nmax = 6;

%%% Autres param�tres
% d�bit de la source [m^3/s]
Q = 1e-3;
% masse volumique du milieu
rho = 1.2;
% c�l�rit� du milieu 
c = 340;

%% Source location

% cartesian coordinates [m]
% xs = [ .3 ; 0 ; -.6];
% ys = [ .3 ; .5 ; 0];
% zs = [ .3 ; -.3 ; -.3];
xs = [ .3 ];
ys = [ .3 ];
zs = [ .3 ];
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
                'ReconstructRadius',a,...
                'c',c,...
                'rho',rho,...
                'Q',Q,...
                'method',sim_method,...
                'incidentOnly',0,...
                'doplot',0);

            
[P] = generateSimu(Rm,Rs,pp_simu);        

%% Affichage 
pp_interp = struct('sphereRadius',a,...
            'numAngle',30,...
            'interpType','natural',...
            'doplot',1);
        
pp_plot = struct('showSource',1,...
                'showMic',1,...
                'plt_typ',plt_typ,...
                'clim',[70 110],...
                'fontSize',12,...
                'freq',f);     
% fig_path = 'C:\Users\Th�o\Documents\1_WORK\01_ENSIM\5A\Projet 5A\Rapports\Fiches de suivi\1912_breve_avancement\';
handle = pressureMeasurementVisu(P,Rm,Rs,pp_interp,pp_plot);
colormap('gray')
view(-45,14);
% fname = sprintf('simu3D_f_%d_view1',f);
% printFigFmt(handle,fig_path,fname,'.png');

handle2 = pressureMeasurementVisu(P,Rm,Rs,pp_interp,pp_plot);
colormap('gray')
view(135,60);
% fname = sprintf('simu3D_f_%d_view2',f);
% printFigFmt(handle2,fig_path,fname,'.png');

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
