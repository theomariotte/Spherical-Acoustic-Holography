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
reg_param = [1e-7 1e-8 1e-9 1e-10 1e-11 1e-12];



% Microphone array radius
a = 15e-2;
% ka product
ka = 0.5;
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

src_loc = [rho_es(:) theta_es(:) phi_es(:)];
src_loc_cart = [X_es(:) Y_es(:) Z_es(:)];

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

%% Compute the strength of each source

pp = struct('c',c,...
            'a',a,...
            'freq',f,...
            'maxOrder',Nmax,...
            'reg_parameter',reg_param);

[q,GN_mat] = getSourceStrength(P_meas,mic_loc,src_loc,pp);
if length(reg_param) == 1
    q_mat = reshape(q,sz);
    p_test = GN_mat * q;
end
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
G = zeros(size(GN_mat));

for ii = 1 : size(q,2)

    % Loop over reconstruction virtual microphones
    for m = 1 : M    
        % location of the current reconstructio mic
        r_m = R_reconstruc(m,:);      
        % Loop over each equivalent source    
        for l = 1 : L
            % Location of the current equivalent source
            r0 = src_loc_cart(l,:);
            R = norm(r_m - r0);
            G(m,l) = exp(-1i * k * R)/(4*pi*R);         
        end
    end

    % incident sound field reconstruction
    p_i = G * q(:,ii);
    
    % test de reconstruction
    pref = 20e-6;
    hh = figure('Name',sprintf('Test de reconstruction mu = %.1e',reg_param(ii)));
    hold on
    plot(flipud(xs - x_tst),flipud(20*log10(abs(p_cmp/pref))),'--k','linewidth',2)
    plot(flipud(xs-x_reconstruc),flipud(20*log10(abs(p_i/pref))),'color',[.5 .5 .5],'linewidth',3);
    hold off
    xlabel('location on the reconstruted axis','interpreter','latex')
    ylabel('SPL [dB]','interpreter','latex')
    grid on
    legend('Simulation','Reconstruction')
    title(['$\mu = $' sprintf('%.1e',reg_param(ii))],'interpreter','latex')
    set(gca,'fontsize',12)
 
    
%     fig_path = 'C:\Users\Théo\Documents\1_WORK\01_ENSIM\5A\Projet 5A\Rapports\Fiches de suivi\1912_breve_avancement\';
%     fname = sprintf('SESM_f_%d',round(f));
%     printFigFmt(hh,fig_path,fname,'eps');
end

%% figures
if length(reg_param) == 1
    % pression simulée sur la sphère (interpolée pour une lecture plus simple)
    pp_interp = struct('sphereRadius',a,...
                'numAngle',30,...
                'interpType','natural',...
                'doplot',0);

    pp_plot = struct('showSource',1,...
                    'showMic',0,...
                    'plt_typ',plt_typ,...
                    'clim',[],...
                    'fontSize',12,...
                    'freq',f);              

    handle1 = pressureMeasurementVisu(P_meas,Rm,Rs,pp_interp,pp_plot);
    colormap('gray')

    hold(gca,'on')
    plot3(X_es,Y_es,Z_es,'o','color',[.5 .5 .5]);
    hold(gca,'off')

    figure('Name','Equivalent sources strength')
    surf(Y_es,Z_es,abs(q_mat))
    colormap('jet')
    shading('interp')
    xlabel('Y')
    ylabel('Z')

    pp_interp = struct('sphereRadius',a,...
                'numAngle',30,...
                'interpType','natural',...
                'doplot',0);

    pp_plot = struct('showSource',1,...
                    'showMic',0,...
                    'plt_typ',plt_typ,...
                    'clim',[],...
                    'fontSize',12,...
                    'freq',f);     

    handle2= pressureMeasurementVisu(p_test,Rm,Rs,pp_interp,pp_plot);
    title('Reconstruction of the sound field on the sphere','interpreter','latex')
end



