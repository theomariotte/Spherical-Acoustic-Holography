% test simulation mesures
clear; clc;
close all;

%% initialisation

% matlabhome = 'H:\Mes documents\5A\Projet_5A\19_10_spherical_NAH\';
matlabhome = 'D:\Users\e152210\19_10_spherical_NAH\';
simu_path = [matlabhome 'simulation\'];
plot_path = [matlabhome 'plot\'];
sph_func_path = [matlabhome 'spherical_functions\'];
data_path = [matlabhome 'data\'];

addpath(plot_path);
addpath(sph_func_path);
addpath(simu_path);
addpath(data_path);

%% Paramètres

% type de calcul : (1) antenne analytique ; (2) vraies positions
typ = 1;

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
% radius of the sphere where the simulation is computed
r_cmp = 2*a;
% fréquence de travail [Hz]
f = 1000;
% Maximum order (Bessel, Hankel, SH)
Nmax = 6;

%%% Autres paramètres
% débit de la source [m^3/s]
Q = 1e-3;
% masse volumique du milieu
rho = 1;
% célérité du milieu 
c = 340;

%% Source location

% cartesian coordinates [m]
xs = 0;
ys = 0;
zs = 0.5;

Rs = [xs ys zs];

%% Load microphone locations
if typ == 1
    
    % grid in azimuth and elevation
    dx = pi/40;
    az_q = (0 : 2*dx : 2*pi)';
    elev_ = (-pi/2 : dx : pi/2)';   
    [az_q,elev_] = meshgrid(az_q,elev_);
    % constant radius
    r_ = r_cmp * ones(size(elev_));
    
    % grid in catesian frame
    [x_,y_,z_] = sph2cart(az_q,elev_,r_);
    sz_ = size(x_);
    
    % microphone's locations
    Rm = [x_(:) y_(:) z_(:)];   
    
    % Interpolation ?
    dointerp = 0;
    
elseif typ == 2
    Nmic = 36;
    filename =  '3Dcam36Officiel.txt';
    fID = fopen([data_path filename],'r');
    sizeA = [Nmic 3];
    mic_loc_tmp = fscanf(fID,'%f');
    fclose(fID);

    Rm = transpose(reshape(mic_loc_tmp,[3 Nmic]));  
    dointerp = 1;
end


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

%% Pressure interpolation to plot on a sphere
if dointerp
    
    % microphone's locations in spherical coordinates        
    x = Rm(:,1);
    y = Rm(:,2);
    z = Rm(:,3);
    
    % interpolation
    F = scatteredInterpolant(x,y,z,P,'natural');

    % grille sur laquelle interpoler
    nTheta0 = 50; 
    nPhi0 = nTheta0;
    az_q = linspace(0, 2*pi, nPhi0);
    elev_q   = linspace(-pi/2, pi/2, nTheta0);
    [az_q, elev_q] = meshgrid(az_q, elev_q);
    
    % rayon de la sphère de calcul
    radius = r_cmp*ones(nTheta0, nPhi0);
    
    % grille d'interpolation en cartésien
    [x_,y_,z_] = sph2cart(az_q, elev_q, radius);
    
    % Pression interpolée sur la nouvelle grille
    P_interp = F(x_,y_,z_);
                  
else    
    P_interp = reshape(P,sz_);
end
%% choose plot
if  strcmp(plt_typ,'real') == 1
    P2plot = real(P_interp);
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
if r_cmp > a
    s1 = surf(x_,y_,z_,P2plot);
    colormap('jet')
    if exist('clim','var')    
        set(gca,'clim',clim);
    end
    shading('interp')
    colorbar
    set(s1,'FaceAlpha',0.5)
    
    [xplt,yplt,zplt] = sphere(30);         
    [az_plt,elev_plt,r_plt] = cart2sph(xplt,yplt,zplt);
    r_plt = a*r_plt;
    [x_plt,y_plt,z_plt] = sph2cart(az_plt,elev_plt,r_plt);
    s1 = surf(x_plt,y_plt,z_plt,ones(size(z_plt)));
    set(s1,'facecolor',[0 0 0]);       
else
    surf(x_,y_,z_,P2plot)    
    colormap('jet')
    if exist('clim','var')    
        set(gca,'clim',clim);
    end
    shading('interp')
    colorbar
end

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
