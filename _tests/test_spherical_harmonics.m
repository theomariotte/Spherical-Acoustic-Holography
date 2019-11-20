%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Script de test du calcul des harmoniques sphériques. Cette grandeur
% intervient lors de la résolution de l'équation d'onde sphérique. Elles
% permettent, dans le cas de l'holographie acoustique sphérique, de
% définir la transformée de Fourier sphérique (projection sur les
% harmoniques sphériques). 
%
% Ce code de test permet de comparer quelques harmoniques sphériques avec
% les résultats théoriques donnés par Williams (paramètre docomp = 1). 
% La fonction de visualisation est également testée. Elle permet de tracer
% les harmoniques sphériques soit en 3D, soit projetées sur une sphère,
% soit projetées sur chaque plan. 
%
% Théo Mariotte
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear variables; clc;
close all;

%% Initialization
set(0, 'defaultAxesTickLabelInterpreter','latex'); 
set(0, 'defaultLegendInterpreter','latex');

% enregistrement images
impath = 'H:\Mes documents\5A\Projet_5A\Images\SphericalHarm\';
ext = 'eps';
doprint = 0;

% type de plot
typ = '3Dplot';

% distance bewteen each angle
dx = pi/200;

% comparaison resultats williams et calcul
docomp = 1;

% Theta : elevation angle
elev = -pi/2 : dx : pi/2;

% Phi : azimuthal angle
az = 0 : dx*2 : 2*pi;

% grid building
[phi,theta] = meshgrid(az,elev);

% Legendre associated function's degree 
n = 2;
m = 2;

% valeur de l'harmonique sphérique pour cette vvaleur d'angle
pp = struct('norm',1,...
            'doplot',1);
[Y,h_legendre] = getSphericalHarmonics(theta,phi,n,m,pp);

if doprint
    fname = sprintf('Legendre_func_n_%02d_m_%02d',n,m);
    flag = printFigFmt(h_legendre,impath,fname,ext);      
end


%% test fonction de visualisation
P = 2;
pp = struct(...
            'pltype',typ,...
            'valtype','imag',...
            'orders',[n m],...
            'fontsize',12);

h_SH = SHvisualization(pp,theta,phi,Y);
colorbar;

pp = struct(...
            'pltype',typ,...
            'valtype','real',...
            'orders',[n m],...
            'fontsize',12);

h_SH = SHvisualization(pp,theta,phi,Y);
colorbar;

if doprint
    fname = sprintf('SH_n_%02d_m_%02d',n,m);
    flag = printFigFmt(h_SH,impath,fname,ext);      
end


%% Quelques harmoniques sphériques théoriques pour tester (voir Williams)

% angle d'élévation corrigé pour être entre 0 et \pi
theta_calc = theta + pi/2;

% cartesian coordinates
if n == 2 && m == 0
    % Y_2^0
    Y_th = sqrt(5/(16*pi)) * (-1 + 3*(cos(theta_calc).^2));
elseif n == 1 && m == -1
    % Y_1^-1
    Y_th = sqrt(3/(8*pi)) * sin(theta_calc).*exp(-1j*phi);
elseif n == 3 && m == 1
    Y_th = sqrt(21/(64*pi)) * exp(1j*phi) .* (1 - 5*(cos(theta_calc)).^2) .* sin(theta_calc);
elseif n == 3 && m == 3
    Y_th = -5/8 * sqrt(7/(5*pi)) * exp(3*1j*phi) .* (sin(theta_calc)).^3;
elseif n == 2 && m == 2
    Y_th = 3 * sqrt(5/(96*pi)) * exp(2*1j*phi) .*(sin(theta_calc)).^2;
else
    docomp = 0;
end

%%% Figure comparison

if docomp
    pp = struct(...
        'pltype',typ,...
        'valtype','imag',...
        'orders',[n m],...
        'fontsize',12);
    h_th = SHvisualization(pp,theta,phi,Y_th);
    title(sprintf('Theoretical spherical harmonic : (%s) $Y_%d^%d$',pp.valtype,n,m),'interpreter','latex')
    colorbar;
    
    pp = struct(...
            'pltype',typ,...
            'valtype','real',...
            'orders',[n m],...
            'fontsize',12);

    h_SH = SHvisualization(pp,theta,phi,Y_th);
    title(sprintf('Theoretical spherical harmonic : (%s) $Y_%d^%d$',pp.valtype,n,m),'interpreter','latex')
    
if doprint
    fname = sprintf('SH_th_n_%02d_m_%02d',n,m);    
    flag = printFigFmt(h_th,impath,fname,ext);    
end
   
colorbar;
end