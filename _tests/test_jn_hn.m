%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Script permettant de tester l'implémentation des fonctions de Bessel et
% Hankel sphériques. La fonction de Hankel sphérique est définie telle que 
%                       hn = jn + i.yn,
% avec jn la fonction de Bessel sphérique du premier type et yn la fonction
% de Bessel sphérique du second type.
% La dérivée de ces fonctions par rapport à la variable r est
% implémantée et testée (elle sera utile pour implémenter S-NAH).
% Le code de visualisation de ces fonctions est également testée.
%
% Théo Mariotte
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear; clc; 
close all;

%% initialisation
matlabhome = 'H:\Mes documents\5A\Projet_5A\19_10_spherical_NAH\';
plot_path = [matlabhome 'plot\'];
sph_func_path = [matlabhome 'spherical_functions\'];

set(0, 'defaultAxesTickLabelInterpreter','latex'); 
set(0, 'defaultLegendInterpreter','latex');

addpath(plot_path);
addpath(sph_func_path);

%% paramètres
dr = .01;
r_in = 0 : dr : 20;
% r_in = 2;
% max function order
n = 6;

% Spherical Bessel function
[hn,dhn,r] = SphericalHankel1(n,r_in);
[jn,djn,~,~,~] = SphericalBessel(n,r_in);

% test validity of the differentiation using Wronksian relation 
% (i.e. jn(r)hn'(r) - jn'(r)hn(r) = i/r^2, i complex)
Gn = jn .* dhn - hn .* djn;

% theoretical value
Gn_th = 1i./(r.^2);

% squared error
e_r = abs(real(Gn_th - Gn)).^2;
e_i = abs(imag(Gn_th - Gn)).^2;

RMSE_r = sqrt(mean(e_r));
RMSE_i = sqrt(mean(e_i));

fprintf('RMSE real part : %.4g\n',RMSE_r)
fprintf('RMSE imaginary part : %.4g\n',RMSE_i)

%% Hankel function (real and imaginary parts) + derivatives

pp = struct('diff_hn',dhn,...
            'order',n,...
            'fontsize',12,...
            'x_label','Radius r');
RadialFuncVisu(hn,r,pp);


%% Comparison Gn and Gn_th

h = figure('Name','Comparison');

subplot(2,1,1)
hold on
plot(r,real(Gn),'k','linewidth',2)
plot(r,real(Gn_th),'k','linewidth',2,'linestyle','--')
plot(r,e_r,'r','linewidth',2)
legend('Estimated','Theoretical','squared error')
grid on
hold off

subplot(2,1,2)
hold on
plot(r,imag(Gn),'k','linewidth',2)
plot(r,imag(Gn_th),'k','linewidth',2,'linestyle','--')
plot(r,e_i,'r','linewidth',2)
set(gca,'ylim',[-200 200])
legend('Estimated','Theoretical','squared error')
grid on
hold off




