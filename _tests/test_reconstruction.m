% test simulation mesures
clear; clc;
close all;

%% Paramètres

data_path = 'data/';
% Méthode de simulation 
% ('Propagator') Utilise le propagateur
% défini par Williams dans 'Intensity vector reconstruction'. 
sim_method = 'Propagator';

%%% Plot parameters
% Tracé de la pression acoustique
% plt_typ : ('real') real part ; ('dB') decibels ; ('abs') absolute value
plt_typ = 'abs';

%%% Paramètres de reconstruction
% fréquence de travail [Hz]
f = 100;
% Maximum order (Bessel, Hankel, SH)
Nmax = 5;
% rayon de reconstruction
r_cmp = 30e-2;


%% Source location
% cartesian coordinates [m]
xs = 0.;
ys = 0.;
zs = 0.5;

Rs = [xs ys zs];

%% Load microphone locations
Nmic = 36;
filename =  '3Dcam36Officiel.txt';
fID = fopen([data_path filename],'r');
sizeA = [Nmic 3];
mic_loc_tmp = fscanf(fID,'%f');
fclose(fID);

Rm = transpose(reshape(mic_loc_tmp,[3 Nmic]));  
dointerp = 0;


%% Génération d'une simulation de mesure de pression sur l'antenne de microphones

%%% paramètre antenne
% Microphone array radius
a = 15e-2;

%%% paramètres du milieu de propagation (air)
% débit de la source [m^3/s]
Q = 1e-3;
% masse volumique du milieu
rho = 1.2;
% célérité du milieu 
c = 340;

pp_simu = struct('MaxOrder',Nmax,...
               'freq',f,...
               'SphereRadius',a,...
               'c',c,...
               'rho',rho,...
               'Q',Q,...
               'method','Propagator',...
               'doplot',0);
            
[P_simu] = generateSimu(Rm,Rs,pp_simu);        

%% Compute the reconstruction on a bigger sphere

[phi_m,theta_m] = cart2sph(Rm(:,1) ,Rm(:,2),Rm(:,3));

pp_reconstruction = struct('freq',f,...
                            'R',r_cmp,...
                            'a',a,...
                            'c',c,...
                            'maxOrder',Nmax,...
                            'incidentOnly',0);

[P_reconstruct] = SNAH(theta_m,phi_m,P_simu,pp_reconstruction);


%% Pressure interpolation to plot on a sphere
if dointerp    
    param = struct('numAngle',50,...
                   'interpType','natural',...
                   'sphereRadius',a);    
    [P_interp,x_grid,y_grid,z_grid] = fieldInterpolation(P_simu,Rm,param);                  
else    
    P_interp = reshape(P_simu,sz_grid);
end

%% choose plot
if  strcmp(plt_typ,'real') == 1
    P2plot = real(P_interp,2);
    min_p = min(min(P_interp));
    max_p = max(max(P_interp));
    abs_lim = max(abs([min_p max_p]));
    clim = [-abs_lim abs_lim];
elseif strcmp(plt_typ,'dB') == 1
    pref = 20e-6;
    P2plot = 20*log10(abs(P_interp)/pref);
    clim = [40 110];
elseif strcmp(plt_typ,'abs') == 1
    P2plot = abs(P_interp);        
end

%% Figures

% pressure field with theoretical sphere
h = figure('Name','Pressure field on the sphere');

hold on
surf(x_grid,y_grid,z_grid,P2plot)    
colormap('jet')
if exist('clim','var')    
    set(gca,'clim',clim);
end
shading('interp')
colorbar

plot3(Rs(:,1),Rs(:,2),Rs(:,3),'ro','linewidth',3)
hold off
xlabel('X','interpreter','latex')
ylabel('Y','interpreter','latex')
zlabel('Z','interpreter','latex')
axis('equal')


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
