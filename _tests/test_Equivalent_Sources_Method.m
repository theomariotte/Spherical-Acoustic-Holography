% test S-ESM method
% clear;
clc;
close all;

%% Paramètres
data_path = 'data/';

% Tracé de la pression acoustique
% plt_typ : ('real') real part ; ('dB') decibels ; ('abs') absolute value
plt_typ = 'dB';
sim_method = 'GN';

% number of equivalent sources
N_src_y = 11;
N_src_z = 11;

% Microphone array radius
a = 15e-2;
% fréquence de travail [Hz]
f = 100;
% Maximum order (Bessel, Hankel, SH)
Nmax = 4;
% Noise amp
nAmp = 0.01;
% color limits
clims = [50 85];

%%% Autres paramètres
% débit de la source [m^3/s]
Q = 5e-4;
% masse volumique du milieu
rho_air = 1.2;
% célérité du milieu 
c = 340;

%% Simulation source location

% cartesian coordinates [m]
xs = 0.5;
ys = 0;
zs = 0;

Rs = [xs ys zs];
[rhos,phis,thetas] = sphericalCoordinates(xs,ys,zs);

%% Equivalent sources locations

% equivalent sources plane location from the simulated source
dx = 5e-2;
dy = 1e-2;
dz = 1e-2;

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
                'c',c,...
                'rho',rho_air,...
                'Q',Q,...
                'method',sim_method,...
                'doplot',0);
          
[P] = generateSimu(Rm,Rs,pp_simu);        

pp_simu = struct('MaxOrder',Nmax,...
                'freq',f,...
                'SphereRadius',a,...
                'c',c,...
                'rho',rho_air,...
                'Q',Q,...
                'method',sim_method,...
                'doplot',0);

[P_cmp] = generateSimu(Rm,Rs,pp_simu);        

P_meas = P + nAmp * randn(size(P));

%% Compute the strength of each source

pp = struct('c',c,...
            'a',a,...
            'freq',f,...
            'maxOrder',Nmax,...
            'reg_parameter',1e-10);

[q,GN_mat] = getSourceStrength(P_meas,mic_loc,src_loc,pp);
q_mat = reshape(q,sz);

p_test = GN_mat * q;

%% Compute the free field Green function for each microphone/source couple
% 
% % grille de reonstruction
% 
% x_grid = -2 : 0.1 : 2;
% y_grid = x_grid;
% [Xg,Yg] = meshgrid(x_grid,y_grid);
% sz = size(Xg);
% Zg = zeros(sz);
% 
% Rm = [Xg(:) Yg(:) Zg(:)];
% 
% Nmic = numel(Xg);
% 
% x = X_es(:);
% y = Y_es(:);
% z = Z_es(:);
% 
% L = N_src_y * N_src_z;
% 
% G = zeros(Nmic,L);
% 
% for l = 1 : L
%     r0 = [x(l) y(l) z(l)];
%     for m = 1 : Nmic
%         r = [Rm(m,1) Rm(m,2) Rm(m,3)];
%         R = norm(r0 - r);
%         
%         G(m,l) = exp(1i * k * R)/(4*pi*R);
%         
%     end    
% end
% 
% % reconstruction du champ incident
% p_i = G * q;
% p_i = reshape(p_i,sz);

%% figures

% pression simulée sur la sphère (interpolée pour une lecture plus simple)
pp_interp = struct('sphereRadius',a,...
            'numAngle',30,...
            'interpType','natural',...
            'doplot',0);
        
pp_plot = struct('showSource',1,...
                'showMic',0,...
                'plt_typ',plt_typ,...
                'clim',clims,...
                'fontSize',12,...
                'freq',f);                 
handle1 = pressureMeasurementVisu(P_meas,Rm,Rs,pp_interp,pp_plot);

hold(gca,'on')
plot3(X_es,Y_es,Z_es,'ko');
hold(gca,'off')

figure('Name','Equivalent sources strength')
surf(Y_es,Z_es,abs(q_mat))
colormap('jet')
xlabel('Y')
ylabel('Z')

pp_interp = struct('sphereRadius',a,...
            'numAngle',30,...
            'interpType','natural',...
            'doplot',0);
        
pp_plot = struct('showSource',1,...
                'showMic',0,...
                'plt_typ',plt_typ,...
                'clim',clims,...
                'fontSize',12,...
                'freq',f);     
            
handle2= pressureMeasurementVisu(p_test,Rm,Rs,pp_interp,pp_plot);

title('Reconstruction of the sound field on the sphere')
% 
% figure('Name','Free field reconstruction')
% % imagesc(x_grid,y_grid,abs(p_i))
% surf(Xg,Yg,abs(p_i))
% xlabel('X')
% ylabel('Y')

