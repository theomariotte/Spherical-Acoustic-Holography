% test S-ESM method
clear; clc;
close all;

%% Paramètres
data_path = 'data/';


% Tracé de la pression acoustique
% plt_typ : ('real') real part ; ('dB') decibels ; ('abs') absolute value
plt_typ = 'dB';
sim_method = 'GN';

% Equivalent sources type (grid or randomly located on a sphere)
ES_typ = 1;

% regularization type ('lsqr : iterative method ; 'tikho' tykhonov reg)
reg_typ = 'lsqr';

% regularization parameter : if 'lsqr' this defines the tolerance (< 1), if
% 'tykho' this defines the so called regularization parameter 
reg_param = 1e-5;
% maximum LSQR iteration number
maxIt = 200;

% Microphone array radius
a = 15e-2;
% ka product
ka = 2;
% Maximum order (Bessel, Hankel, SH)
Nmax = 4;
% Noise amp
nAmp = 0.001;
% color limits
clims = [40 70];

%%% Autres paramètres

% débit de la source [m^3/s]
Q = 1e-5;
% masse volumique du milieu
rho_air = 1.2;
% célérité du milieu 
c = 340;
% fréquence de travail [Hz]
f = ka*c/(2*pi*a);
%f = 1.1796e+03;
w = 2*pi*f;
k = w/c;
phasefl = 0;

%% Simulation source location

% cartesian coordinates [m]
xs = [0.2];
ys = [0.0];
zs = [0.];

Rs = [xs ys zs];
[rhos,phis,thetas] = sphericalCoordinates(xs,ys,zs);

%% Equivalent sources locations

if ES_typ

    N_src_y = 41;
    N_src_z = 41;
    
    % equivalent sources plane location from the simulated source
    dx = .5e-2;
    dy = 0.1e-2;
    dz = 0.1e-2;

    y_loc_inf = -floor(N_src_y/2) * dy;
    y_loc_sup =  floor(N_src_y/2) * dy;
    z_loc_inf = -floor(N_src_z/2) * dz;
    z_loc_sup =  floor(N_src_z/2) * dz;

    Y_es_tmp = y_loc_inf : dy : y_loc_sup;
    Z_es_tmp = z_loc_inf : dz : z_loc_sup;

    [Y_es,Z_es] = meshgrid(Y_es_tmp,Z_es_tmp);

    X_es = (xs(1)+dx) * ones(size(Y_es));
    sz = size(Y_es);

    [rho_es,phi_es,theta_es] = sphericalCoordinates(X_es,Y_es,Z_es);

    ES_loc = [rho_es(:) theta_es(:) phi_es(:)];
    ES_loc_cart = [X_es(:) Y_es(:) Z_es(:)];

else
    % DOES NOT WORK 
    Ns = 10;
    R = xs + 1e-2;
    [Xs,Ys,Zs] = randomESLocation(R,Ns,1);
    [X_es,Y_es,Z_es] = meshgrid(Xs,Ys,Zs);
    [rho_es,phi_es,theta_es] = sphericalCoordinates(X_es,Y_es,Z_es);
    ES_loc = [rho_es(:) theta_es(:) phi_es(:)];
    ES_loc_cart = [X_es(:) Y_es(:) Z_es(:)];
    
end

%% Microphones location in cartesian frame

Nmic = 36;
filename =  '3Dcam36Officiel.txt';
fID = fopen([data_path filename],'r');
sizeA = [Nmic 3];
mic_loc_tmp = fscanf(fID,'%f');
fclose(fID);

Rm = transpose(reshape(mic_loc_tmp,[3 Nmic]));  

[rho_mic,phi_mic,theta_mic] = sphericalCoordinates(Rm(:,1),Rm(:,2),Rm(:,3));

mic_loc = [rho_mic theta_mic phi_mic];


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
                'doplot',1);
          
[P] = generateSimu(Rm,Rs,pp_simu,phasefl);        

P_meas = P + nAmp * randn(size(P));

%% Compute the transfert matrix based on the Green Neumann function 
Gn = getGreenNeumannFRF(ES_loc,mic_loc,Nmax,k,a);

%% Solve the inverse problem using regularization

% solve ill conditionned problem
q = solveIllPosedProblem(Gn,P_meas,reg_typ,reg_param);    

%% Targeted sound pressure on the reconstruciton axis

% grille de reconstruction sur l'axe radial
r_start = 0.01;
da = 0.001;
x_tst =(r_start : da : xs(1))' ;
x_analytic = Rs(1)-x_tst;
y_analytic = Rs(2)-zeros(size(x_analytic));
z_analytic =  Rs(3)-zeros(size(x_analytic));
rm = [x_analytic y_analytic z_analytic];
R_cmp = sqrt(sum(rm.^2,2));
fact = 1i*w*rho_air*Q;
p_cmp = fact * exp(-1i*k*R_cmp)./(4*pi*R_cmp);

%% Compute the free field Green function for each microphone/source couple

% Reconstruction grid on the radial axis
dr = 0.001;
x_reconstruc = (r_start : dr : xs(1))';
y_reconstruc = zeros(size(x_reconstruc));
z_reconstruc = zeros(size(x_reconstruc));

R_reconstruc = [x_reconstruc y_reconstruc z_reconstruc];

G = getFreeFieldFRF(ES_loc_cart,R_reconstruc,k);

% incident sound field reconstruction
p_i = G * q;

%% reconstruction on a plane 

dx = 1E-2;
dz = 1E-2;

xMin = -a;
xMax = xs;
zMin = -2*a;
zMax = 2*a;

[fieldPointsX,fieldPointsZ] = meshgrid(xMin:dx:xMax,zMin:dz:zMax);
sZ = size(fieldPointsX);
fieldPointsY = zeros(sZ);

fieldPoints = [fieldPointsX(:) fieldPointsY(:) fieldPointsZ(:)];
[rFP,phiFP,thetaFP] = sphericalCoordinates(fieldPointsX(:),fieldPointsY(:),fieldPointsZ(:));
fieldPointsSph = [rFP thetaFP phiFP];

FRF_reconstruct = getFreeFieldFRF(ES_loc_cart,fieldPoints,k);

%FRF_reconstruct = getGreenNeumannFRF(ES_loc,fieldPointsSph,Nmax,k,a);

p_i_far = FRF_reconstruct*q;

p_i_far = reshape(p_i_far,sZ);


%% figures

% test de reconstruction
pref = 20e-6;
if length(xs) < 2    
    hh = figure('Name','Test de reconstruction LSQR');
    hold on
    plot(flipud(xs - x_tst),flipud(20*log10(abs(p_cmp/pref))),'--k','linewidth',2)
    plot(flipud(xs-x_reconstruc),flipud(20*log10(abs(p_i/pref))),'color',[.5 .5 .5],'linewidth',3);
    hold off
    xlabel('Distance to the source [m]','interpreter','latex')
    ylabel('SPL [dB]','interpreter','latex')
    grid on
    legend('Simulation','Reconstruction')
    title('Reconstruction using LSQR')
    set(gca,'fontsize',12)
end

if ES_typ
    
    % for testing
    q_mat = reshape(q,sz);

    figure('Name','Equivalent sources strength')
    imagesc(Y_es_tmp,Z_es_tmp,abs(q_mat),'Interpolation','bilinear')
    %shading  interp
    colormap jet
    xlabel('Y')
    ylabel('Z')

end

p_test = Gn * q;
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
title('Reconstruction of the sound field on the sphere','interpreter','latex')
    
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

handle3= pressureMeasurementVisu(P_meas,Rm,Rs,pp_interp,pp_plot);
title('Sound pressure measurement (interpolated)','interpreter','latex')


figure('Name','Field points');
hold on
imagesc(xMin:dx:xMax,zMin:dz:zMax,20*log10(abs(p_i_far/pref)),'Interpolation','bilinear')
contour(xMin:dx:xMax,zMin:dz:zMax,20*log10(abs(p_i_far/pref)),10,'linecolor',[0 0 0])
ll = rectangle('position',[-a -a 2*a 2*a],'Curvature',[1 1]);
set(ll,'linewidth',3)
hold off
xlabel('X')
ylabel('Z')
title('Reconstruction of the incident pressure over a plane grid')
axis equal

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





