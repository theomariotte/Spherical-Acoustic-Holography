% test S-ESM method
clear;
clc;
close all;

%% Paramètres
data_path = 'data/';

% Tracé de la pression acoustique
% plt_typ : ('real') real part ; ('dB') decibels ; ('abs') absolute value
plt_typ = 'dB';
sim_method = 'GN';

% number of equivalent sources
N_src_y = 15;
N_src_z = 15;

% regularization parameter
num_reg = 500;
reg_param = linspace(1e-5,10,num_reg);

% Microphone array radius
a = 15e-2;
% ka product
ka = 2;
% Maximum order (Bessel, Hankel, SH)
Nmax = 4;
% Noise amp
nAmp = 0.0;
% color limits
clims = [50 100];

%%% Autres paramètres

% débit de la source [m^3/s]
Q = 5e-4;
% masse volumique du milieu
rho_air = 1.2;
% célérité du milieu 
c = 340;
% fréquence de travail [Hz]
f = ka*c/(2*pi*a);
w = 2*pi*f;
k = w/c;

%% Simulation source location

% cartesian coordinates [m]
xs = 0.8;
ys = 0;
zs = 0;

Rs = [xs ys zs];
[rhos,phis,thetas] = sphericalCoordinates(xs,ys,zs);

%% Equivalent sources locations

% equivalent sources plane location from the simulated source
dx = 1e-2;
dy = 0.5e-2;
dz = 0.5e-2;

y_loc_inf = -floor(N_src_y/2) * dy;
y_loc_sup =  floor(N_src_y/2) * dy;
z_loc_inf = -floor(N_src_z/2) * dz;
z_loc_sup =  floor(N_src_z/2) * dz;

Y_es_tmp = y_loc_inf : dy : y_loc_sup;
Z_es_tmp = z_loc_inf : dz : z_loc_sup;

[Y_es,Z_es] = meshgrid(Y_es_tmp,Z_es_tmp);

X_es = (xs+dx) * ones(size(Y_es));
sz = size(Y_es);

[rho_es,phi_es,theta_es] = sphericalCoordinates(X_es,Y_es,Z_es);

ES_loc = [rho_es(:) theta_es(:) phi_es(:)];
ES_loc_cart = [X_es(:) Y_es(:) Z_es(:)];

%% Chargement ds coordonnées des microphones (cartésien)

Nmic = 36;
filename =  '3Dcam36Officiel.txt';
fID = fopen([data_path filename],'r');
sizeA = [Nmic 3];
mic_loc_tmp = fscanf(fID,'%f');
fclose(fID);

Rm = transpose(reshape(mic_loc_tmp,[3 Nmic]));  

[rho,phi,theta] = sphericalCoordinates(Rm(:,1),Rm(:,2),Rm(:,3));

mic_loc = [rho theta phi];


%% calculate the soundfield at each microphone

pp_simu = struct('MaxOrder',Nmax,...
                'freq',f,...
                'SphereRadius',a,...
                'ReconstructRadius',a,...
                'c',c,...
                'rho',rho_air,...
                'Q',Q,...
                'method',sim_method,...
                'incidentOnly',0,...
                'doplot',0);
          
[P] = generateSimu(Rm,Rs,pp_simu);        

P_meas = P + nAmp * randn(size(P));

%% Compute the transfert matrix based on the Green Neumann function 
Gn = getGreenNeumannFRF(ES_loc,mic_loc,Nmax,k,a);

%% Solve the inverse problem using regularization
GCV = zeros(num_reg,1);
for i_lambda = 1 : num_reg
    
    % solve ill conditionned problem
    lambda = reg_param(i_lambda);
    q = solveIllPosedProblem(Gn,P_meas,lambda);

    % for testing
    q_mat = reshape(q,sz);
    p_test = Gn * q;

    %% Generalised cross validation to get the optimal regularisation parameter
    lambda = reg_param(i_lambda);
    G_prod = Gn'*Gn ;
    I = eye(size(G_prod));
    inv_mat = inv(G_prod + lambda * I);
    B = Gn * (inv_mat * Gn');
    den = (1/Nmic * sum( diag( eye(size(B)) -B ) ) )^2;
    GCV(i_lambda) = 1/Nmic * norm(1 - B*P_meas)^2 / den;
    
end


h = figure('Name','GCV');
plot(reg_param,abs(GCV),'k','linewidth',2);
xlabel('Regularization parameter \mu');
ylabel('GCV')
grid on


return;
%% Targeted sound pressure on the reconstrucito axis

% grille de reconstruction sur l'axe radial
r_start = 0.;
da = 0.001;
x_tst =(r_start : da : xs - da)' ;
x_analytic = x_tst - Rs(1);
y_analytic = zeros(size(x_analytic)) - Rs(2);
z_analytic = zeros(size(x_analytic)) - Rs(3);
rm = [x_analytic y_analytic z_analytic];
R_cmp = sqrt(sum(rm.^2,2));
G_cmp = exp(-1i*k*R_cmp)./(4*pi*R_cmp);
p_cmp = 1i*(2*pi*f)*rho_air*Q*G_cmp;

%% Compute the free field Green function for each microphone/source couple

% Grille de reconstruction sur l'axe radial
dr = 0.01;
x_reconstruc = (r_start : dr : xs)';
y_reconstruc = zeros(size(x_reconstruc));
z_reconstruc = zeros(size(x_reconstruc));

R_reconstruc = [x_reconstruc y_reconstruc z_reconstruc];

M = length(x_reconstruc);
L = size(q,1);
G = zeros(size(Gn));

% Loop over reconstruction virtual microphones
for m = 1 : M    
    % location of the current reconstructio mic
    r_m = R_reconstruc(m,:);      
    % Loop over each equivalent source    
    for l = 1 : L
        % Location of the current equivalent source
        r0 = ES_loc_cart(l,:);
        R = norm(r_m - r0);
        G(m,l) = exp(-1i * k * R)/(4*pi*R);         
    end
end

% incident sound field reconstruction
p_i = G * q;

% test de reconstruction
pref = 20e-6;
hh = figure('Name',sprintf('Test de reconstruction mu = %.1e',reg_param));
hold on
plot(flipud(xs - x_tst),flipud(20*log10(abs(p_cmp/pref))),'--k','linewidth',2)
plot(flipud(xs-x_reconstruc),flipud(20*log10(abs(p_i/pref))),'color',[.5 .5 .5],'linewidth',3);
hold off
xlabel('Distance to the source [m]','interpreter','latex')
ylabel('SPL [dB]','interpreter','latex')
grid on
legend('Simulation','Reconstruction')
title(['$\mu = $' sprintf('%.1e',reg_param)],'interpreter','latex')
set(gca,'fontsize',12)
 
    
%     fig_path = 'C:\Users\Théo\Documents\1_WORK\01_ENSIM\5A\Projet 5A\Rapports\Fiches de suivi\1912_breve_avancement\';
%     fname = sprintf('SESM_f_%d',round(f));
%     printFigFmt(hh,fig_path,fname,'eps');





