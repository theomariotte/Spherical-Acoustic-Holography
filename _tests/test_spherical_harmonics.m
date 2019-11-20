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

% valeur de l'harmonique sph�rique pour cette vvaleur d'angle
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


%% Quelques harmoniques sph�riques th�oriques pour tester (voir Williams)

% angle d'�l�vation corrig� pour �tre entre 0 et \pi
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