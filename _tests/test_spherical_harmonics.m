%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Script de test du calcul des harmoniques sph�riques. Cette grandeur
% intervient lors de la r�solution de l'�quation d'onde sph�rique. Elles
% permettent, dans le cas de l'holographie acoustique sph�rique, de
% d�finir la transform�e de Fourier sph�rique (projection sur les
% harmoniques sph�riques). 
%
% Ce code de test permet de comparer quelques harmoniques sph�riques avec
% les r�sultats th�oriques donn�s par Williams (param�tre docomp = 1). 
% La fonction de visualisation est �galement test�e. Elle permet de tracer
% les harmoniques sph�riques soit en 3D, soit projet�es sur une sph�re,
% soit projet�es sur chaque plan. 
%
% Th�o Mariotte
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear variables; clc;
close all;

%% Initialization

matlabhome = 'H:\Mes documents\5A\Projet_5A\19_10_spherical_NAH\';
plot_path = [matlabhome 'plot\'];
sph_func_path = [matlabhome 'spherical_functions\'];

set(0, 'defaultAxesTickLabelInterpreter','latex'); 
set(0, 'defaultLegendInterpreter','latex');

addpath(plot_path);
addpath(sph_func_path);

% distance bewteen each angle
dx = pi/200;

% comparaison resultats williams et calcul
docomp = 1;

% Theta : elevation angle
elev = 0 : dx : pi;
% elev = -pi/2 : dx : pi/2;
% elev = -pi : dx : pi;
% Phi : azimuthal angle
az = 0 : dx*2 : 2*pi;
% grid building
[phi,theta] = meshgrid(az,elev);

% Legendre associated function's degree 
n = 4;
m = 1;

% valeur de l'harmonique sph�rique pour cette vvaleur d'angle
pp = struct('norm',1,...
            'doplot',0);
[Y] = getSphericalHarmonics(theta,phi,n,m,pp);

%% test fonction de visualisation
P = 2;
pp = struct(...
            'pltype','2Dplot',...
            'valtype','real',...
            'orders',[n m],...
            'fontsize',12);

h_SH = SHvisualization(pp,theta,phi,Y);



%% Quelques harmoniques sph�riques th�oriques pour tester (voir Williams)
% cartesian coordinates

if n == 2 && m == 0
    % Y_2^0
    Y_th = sqrt(5/(16*pi)) * (-1 + 3*(cos(theta).^2));
elseif n == 1 && m == -1
    % Y_1^-1
    Y_th = sqrt(3/(8*pi)) * sin(theta).*exp(-1j*phi);
elseif n == 3 && m == 1
    Y_th = sqrt(21/(64*pi)) * exp(1j*phi) .* (1 - 5*(cos(theta)).^2) .* sin(theta);
elseif n == 3 && m == 3
    Y_th = -5/8 * sqrt(7/(5*pi)) * exp(3*1j*phi) .* (sin(theta)).^3;
else
    docomp = 0;
end

%%% Figure comparison

if docomp
    
    h2 = figure('Name','Theoretical SH');
    
    % cartesian coordinates of the theoretical SH
    [Xth,Yth,Zth] = sph2cart(phi,theta-pi/2,real(Y_th));
    
    surf(Xth,Yth,Zth,real(Y_th))
    shading('interp')
    title(sprintf('Theoretical spherical harmonic $Y_{%d}^{%d}$',n,m),'interpreter','latex')
%     set(gca,'xlim',[-0.5 0.5],'ylim',[-0.5 0.5],'zlim',[-0.5 0.5])
    colormap('jet')
    axis equal
end