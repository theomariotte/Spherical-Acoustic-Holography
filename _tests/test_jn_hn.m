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
set(0, 'defaultAxesTickLabelInterpreter','latex'); 
set(0, 'defaultLegendInterpreter','latex');
doprint = 1;
impath = 'H:\Mes documents\5A\Projet_5A\Images\';
%% paramètres
dr = .01;
r_in = 0 : dr : 20;
% r_in = 2;
% max function order
n = 5;

% Spherical Bessel function
[hn,dhn,r] = SphericalHankel1(n,r_in);
[jn,djn,yn,dyn,~] = SphericalBessel(n,r_in);

% test validity of the differentiation using Wronksian relation 
% (i.e. jn(r)hn'(r) - jn'(r)hn(r) = i/r^2, i complex)
Gn = jn .* dhn - hn .* djn;
Gn_2 = jn .* dyn - yn .* djn;

% theoretical value
Gn_th = 1i./(r(1:2:end).^2);
Gn_th2 = 1./(r(1:2:end).^2);

% squared error
% e_r = abs(Gn_th2 - Gn_2).^2;
% e_i = abs(imag(Gn_th - Gn)).^2;

% RMSE_r = sqrt(mean(e_r));
% RMSE_i = sqrt(mean(e_i));

% fprintf('RMSE real part : %.4g\n',RMSE_r)
% fprintf('RMSE imaginary part : %.4g\n',RMSE_i)

%% Hankel function (real and imaginary parts) + derivatives

pp = struct('diff_hn',dhn,...
            'order',n,...
            'fontsize',12,...
            'x_label','Radius r');
h1 = RadialFuncVisu(hn,r,pp);

if doprint
   fname = sprintf('Bessel_order_%02d',n);
    outName = getOutFileName(impath,fname,'.eps');
    print(h1,outName,'-depsc'); 
end

%% Comparison Gn and Gn_th

h = figure('Name','Bessel 1 et 2');

hold on
plot(r(1:2:end),Gn_th2,'k','linewidth',2,'markersize',7);
plot(r,Gn_2,'rx','linewidth',1);
hold off

set(gca,'yscale','log','ylim',[1e-1 1e3],'xlim',[0 3])
set(gca,'fontsize',12);
grid on

title('Wronskien : $j_n(r)$ et $y_n(r)$','interpreter','latex')
xlabel('Rayon $r$','interpreter','latex');
ylabel("$G_{n}(r) = j_n(r)y_n'(r) - j_n'(r)y_n(r)$",'interpreter','latex')
legend('Theorie','Calcul')


if doprint    
    fname = 'Wronskien_jn_yn';

    outName = getOutFileName(impath,fname,'.eps');
    print(h,outName,'-depsc');
end
%%%

h4 = figure('Name','Bessel Hankel');

hold on
plot(r(1:2:end),imag(Gn_th),'k','linewidth',2,'markersize',7)
plot(r,imag(Gn),'rx','linewidth',1)
hold off

set(gca,'yscale','log')
set(gca,'fontsize',12);
set(gca,'ylim',[1e-1 1000],'xlim',[0 3])
grid on

title('Wronskien : $j_n(r)$ et $h_n(r)$','interpreter','latex')
xlabel('Rayon $r$','interpreter','latex');
ylabel("$G_{n}(r) = j_n(r)h_n'(r) - j_n'(r)h_n(r)$",'interpreter','latex')
legend('Theorie','Calcul')

if doprint
    fname = 'Wronskien_jn_hn';

    outName = getOutFileName(impath,fname,'.eps');
    print(h4,outName,'-depsc');
end




