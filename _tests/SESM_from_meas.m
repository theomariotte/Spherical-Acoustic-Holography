% Reconstruction test from experimental data

clear; clc; 
close all;

%% Parameters

fig_path = '/Users/theomariotte/Documents/01_work/ENSIM/5A/Projet/Figures/';

% frequency 
f0 = 500;
c = 340;
H = 1.116;
% marameters for equivalent sources

% location of the source to be characterized
ys = 41.5e-2;
dy = .5e-2;

% distance between equivalent sources on the plabe X,Z
dx = 0.01;
dz = 0.01;

% Taking ground into account (0) no, (1) yes
groundfl = 0;

% location of the source 

% % first case : 7th location
Rs = [0 ys 0];

% second case : 1st and 12th location
%Rs = [[-0.25 ys 0.253];[0.25 ys 0]];

%% load data measured by the array

% pressure measured by the spherical array with 2 sources
dataPath = ['/Users/theomariotte/Documents/01_work/ENSIM/5A/Projet/'...
    'Spherical-Acoustic-Holography/data/12022020_110356/'];

fname = 'data/12022020_110356/Pression_acoustique.mat';

% %pressure measured by the spherical array with 2 sources
% dataPath = ['/Users/theomariotte/Documents/01_work/ENSIM/5A/Projet/'...
%     'Spherical-Acoustic-Holography/data/12022020_111032/'];
% 
% fname = 'data/12022020_111032/Pression_acoustique.mat';

type_data=0; % Pressions acoustiques

no_devices=1:3; % n°carte
no_voies{1}=1:16; % n°voie
no_voies{2}=1:16;
no_voies{3}=1:4;
nbvoies=36;

addpath(dataPath);

[pArrayTime,t,fs]=convert_data_TDMS(nbvoies,type_data,no_devices,no_voies,0);

% reference sound pressure
dataPath = ['/Users/theomariotte/Documents/01_work/ENSIM/5A/Projet/'...
    'Spherical-Acoustic-Holography/data/12022020_112037/'];

addpath(dataPath);

no_devices=3; % n°carte
no_voies{1}=5; % n°voie
nbvoies2=1;

[pRefTime,tref,~]=convert_data_TDMS(nbvoies2,type_data,no_devices,no_voies,0);

% background noise


%% compute spectrum
cd(['/Users/theomariotte/Documents/01_work/'...
    'ENSIM/5A/Projet/Spherical-Acoustic-Holography/']);
M = 2^13;
pref = 20e-6;

Pxy = zeros(M/2+1,nbvoies);
pMeas = zeros(nbvoies,1);

ff = linspace(0,fs/2,M/2+1);
[~,idxF] = min(abs(ff-f0));

% figure('Name','PSD')
% hold on

for k = 1 : nbvoies
   
    crnt_x = eval(sprintf('pArrayTime.voie%d',k));
    Pxy(:,k) = cpsd(pRefTime.voie1,crnt_x,hann(M),M/2,M);
    Pxx = pwelch(pRefTime.voie1,hann(M),M/2,M);
    pMeas(k) = Pxy(idxF,k);%/(sqrt(abs(Pxx(k))^2));
    %plot(ff,10*log10(abs(Pxy(:,k)/(pref^2))));
    
end

handle3 = figure('Name','PSD');
hold on
plot(ff,10*log10(abs(Pxy(:,1)/pref^2)),'color',[.2 .4 .8])
plot(ff,10*log10(abs(Pxy(:,25)/pref^2)),'color',[.8 .3 .1])
hold off
legend('Mic 1','Mic 25');
xlabel('frequency [Hz]')
ylabel('Amplitude [dB]')
grid on
set(gca,'fontsize',12,'xlim',[0 fs/2]);

% fname = 'CPSDexample';
% printFigFmt(handle3,fig_path,fname,'eps');

% grid on
% set(gca,'xlim',[0 fs/2]);

%% compute equivalent sources locations


x_loc_inf = -0.5;
x_loc_sup =  0.5;
z_loc_inf = -0.5;
z_loc_sup =  0.5;

Y_es_tmp = x_loc_inf : dx : x_loc_sup;
Z_es_tmp = z_loc_inf : dz : z_loc_sup;

[X_es,Z_es] = meshgrid(Y_es_tmp,Z_es_tmp);

Y_es = (ys(1)+dy) * ones(size(X_es));
sz = size(X_es);

[rho_es,phi_es,theta_es] = sphericalCoordinates(X_es,Y_es,Z_es);

ES_loc = [rho_es(:) theta_es(:) phi_es(:)];
ES_loc_cart = [X_es(:) Y_es(:) Z_es(:)];

%% Microphones location in cartesian frame
dataPath = ['/Users/theomariotte/Documents/01_work/ENSIM/5A/Projet/'...
    'Spherical-Acoustic-Holography/data/'];

Nmic = 36;
filename =  '3Dcam36Officiel.txt';
fID = fopen([dataPath filename],'r');
sizeA = [Nmic 3];
mic_loc_tmp = fscanf(fID,'%f');
fclose(fID);

Rm = transpose(reshape(mic_loc_tmp,[3 Nmic]));  

% array center was at 1,12M from the ground
Rm(:,3) = Rm(:,3);

[rho_mic,phi_mic,theta_mic] = sphericalCoordinates(Rm(:,1),Rm(:,2),Rm(:,3));

mic_loc = [rho_mic theta_mic phi_mic];


%% compute sources strenght

% array parameters
a = 15e-2;
Nmax = 4;
k = 2*pi*ff(idxF)/c;

pp_interp = struct('sphereRadius',a,...
            'numAngle',30,...
            'interpType','natural',...
            'doplot',0);

pp_plot = struct('showSource',1,...
                'showMic',0,...
                'plt_typ','dB',...
                'clim',[],...
                'fontSize',12,...
                'freq',f0);              

handle1 = pressureMeasurementVisu(pMeas,Rm,Rs,pp_interp,pp_plot);
colormap('jet')
xlabel('Y')
ylabel('X')
zlabel('Z')

%%% transfert function between sources and array 

if groundfl
    
    % Image equivalent sources
    IES_loc_cart = ES_loc_cart;
    IES_loc_cart(:,3) = IES_loc_cart(:,3) - 2*H;
    
    IES_loc_cart = flipud(IES_loc_cart);
    
    [rhoTmp,phiTmp,thetaTmp] = sphericalCoordinates(IES_loc_cart(:,1),...
        IES_loc_cart(:,2),IES_loc_cart(:,3));
    
    IES_loc = [rhoTmp,thetaTmp,phiTmp];
    
    GN_es = getGreenNeumannFRF(ES_loc,mic_loc,Nmax,k,a);
    GN_image = getGreenNeumannFRF(IES_loc,mic_loc,Nmax,k,a);
    
    GN = GN_es + GN_image;
    
else
    
    % green neumann without taking the ground into account
    GN = getGreenNeumannFRF(ES_loc,mic_loc,Nmax,k,a);
    
end

q = solveIllPosedProblem(GN,pMeas,'tikhonov',1e-2);

% % L curve
% pp = 2000;
% lambdaVec = linspace(1e-12,1e-4,pp);
% for k = 1 : pp
%     q = solveIllPosedProblem(GN,pMeas,'tikhonov',lambdaVec(k));
%     normQ(k) = norm(q,2);
%     normDiff(k) = norm(GN*q-pMeas,2);
% end
% 
% figure
% plot(log10(normDiff),log10(normQ),'k')
% xlabel('p-Gq')
% ylabel('q')
% grid on
% 
% return


%% reconstruction on da plane YZ

dx = 1E-2;
dz = 1E-2;

yMin = -a;
yMax = ys;
zMin = -H;
zMax = 0.5;

[fieldPointsY,fieldPointsZ] = meshgrid(yMin:dx:yMax,zMin:dz:zMax);
sZ = size(fieldPointsY);
fieldPointsX = zeros(sZ);

fieldPoints = [fieldPointsX(:) fieldPointsY(:) fieldPointsZ(:)];
[rFP,phiFP,thetaFP] = sphericalCoordinates(fieldPointsX(:),...
    fieldPointsY(:),fieldPointsZ(:));
fieldPointsSph = [rFP thetaFP phiFP];

FRF_reconstruct = getFreeFieldFRF(ES_loc_cart,fieldPoints,k);

%FRF_reconstruct = getGreenNeumannFRF(ES_loc,fieldPointsSph,Nmax,k,a);

p_i_farYZ = FRF_reconstruct*q;

p_i_farYZ = reshape(p_i_farYZ,sZ);

%% reconstruction on da plane XY

dx = 1E-2;
dy = 1E-2;

yMin = -a;
yMax = ys;
xMin = -0.5;
xMax = 0.5;

[fieldPointsX,fieldPointsY] = meshgrid(xMin:dx:xMax,yMin:dy:yMax);
sZ = size(fieldPointsX);
fieldPointsZ = zeros(sZ);

fieldPoints = [fieldPointsX(:) fieldPointsY(:) fieldPointsZ(:)];
[rFP,phiFP,thetaFP] = sphericalCoordinates(fieldPointsX(:),...
    fieldPointsY(:),fieldPointsZ(:));
fieldPointsSph = [rFP thetaFP phiFP];

FRF_reconstruct = getFreeFieldFRF(ES_loc_cart,fieldPoints,k);

%FRF_reconstruct = getGreenNeumannFRF(ES_loc,fieldPointsSph,Nmax,k,a);

p_i_farXY = FRF_reconstruct*q;

p_i_farXY = reshape(p_i_farXY,sZ);

%% reconstruction on a grid in front of the sources plane

% dx = 1E-2;
% dz = 1E-2;
% 
% [fieldPointsX,fieldPointsZ] = meshgrid(xMin:dx:xMax,zMin:dz:zMax);
% sZ = size(fieldPointsX);
% fieldPointsY = zeros(sZ);
% 
% fieldPoints = [fieldPointsX(:) fieldPointsY(:) fieldPointsZ(:)];
% [rFP,phiFP,thetaFP] = sphericalCoordinates(fieldPointsX(:),...
%     fieldPointsY(:),fieldPointsZ(:));
% fieldPointsSph = [rFP thetaFP phiFP];
% 
% FRF_reconstruct = getFreeFieldFRF(ES_loc_cart,fieldPoints,k);
% 
% %FRF_reconstruct = getGreenNeumannFRF(ES_loc,fieldPointsSph,Nmax,k,a);
% 
% p_i_farXZ = FRF_reconstruct*q;
% 
% p_i_farXZ = reshape(p_i_farXZ,sZ);


%% Figures

pp_interp = struct('sphereRadius',a,...
            'numAngle',30,...
            'interpType','natural',...
            'doplot',0);

pp_plot = struct('showSource',0,...
                'showMic',0,...
                'plt_typ','dB',...
                'clim',[],...
                'fontSize',12,...
                'freq',f0);              

handle1 = pressureMeasurementVisu(pMeas,Rm,Rs,pp_interp,pp_plot);
colormap('jet')
xlabel('Y')
ylabel('X')
zlabel('Z')
hold(gca,'on')
plot3(X_es,Y_es,Z_es,'o','color',[.5 .5 .5]);
if groundfl
   plot3(IES_loc_cart(:,1),IES_loc_cart(:,2),IES_loc_cart(:,3),...
       'o','color',[.5 .5 .5]); 
end
[xGround,yGround] = meshgrid(-1:.5:1);
zGround = zeros(size(xGround)) - H;
surf(xGround,yGround,zGround);
hold(gca,'off')
set(gca,'fontsize',20)
view(142,20); 
% 
% fname = sprintf('meas3D_f_%d_view1',f0);
% printFigFmt(handle1,fig_path,fname,'.eps');

%%%%%


q_mat = reshape(q,sz);

hu = figure('Name','Equivalent sources strength');
hold on
imagesc(Y_es_tmp,Z_es_tmp,abs(q_mat),'Interpolation','bilinear')
set(gca,'yDir','normal')
dessin_13HP(0,0.05)
set(gca,'xlim',[min(Y_es_tmp) max(Y_es_tmp)],'ylim',[min(Z_es_tmp) max(Z_es_tmp)])
hold off
colormap jet
xlabel('Y')
ylabel('Z')
title('Force des sources équivalentes')
set(gca,'fontsize',20)

% fname = sprintf('sourcesStrength');
% printFigFmt(hu,fig_path,fname,'.eps');

%%%%

h1 = figure('Name','Field points');
hold on
imagesc(yMin:dx:yMax,zMin:dz:zMax,20*log10(abs(p_i_farYZ/pref)),...
    'Interpolation','bilinear')
contour(yMin:dx:yMax,zMin:dz:zMax,20*log10(abs(p_i_farYZ/pref)),10,...
    'linecolor',[0 0 0])
ll = rectangle('position',[-a -a 2*a 2*a],'Curvature',[1 1]);
set(ll,'linewidth',3)
l1 = line([yMax yMax],[-0.375 0.375],'linewidth',3,'color',...
    [0 0 0],'linestyle','-');
hold off
xlabel('X')
ylabel('Z')
%title('Reconstruction du champ dans le plan XZ')
axis equal
cc = colorbar;
cc.Label.String = 'dB SPL';

colormap('jet')
set(gca,'clim',[40 100],'fontsize',20)
fname = sprintf('MeasReconstructionPlaneYzGround%d',groundfl);
printFigFmt(h1,fig_path,fname,'eps');

%%%%

h2 = figure('Name','Field points');
hold on
imagesc(xMin:dx:xMax,yMin:dy:yMax,20*log10(abs(p_i_farXY/pref)),...
    'Interpolation','bilinear')
contour(xMin:dx:xMax,yMin:dy:yMax,20*log10(abs(p_i_farXY/pref)),...
    10,'linecolor',[0 0 0])
ll = rectangle('position',[-a -a 2*a 2*a],'Curvature',[1 1]);
set(ll,'linewidth',3)
l1 = line([-0.375 0.375],[yMax yMax],'linewidth',3,...
    'color',[0 0 0],'linestyle','-');
hold off
xlabel('Y')
ylabel('X')
%title('Reconstruction du champ dans le plan XY')
axis equal
cc = colorbar;
cc.Label.String = 'dB SPL';
colormap('jet')
set(gca,'clim',[40 100],'fontsize',20)

fname = sprintf('MeasReconstructionPlaneXYGround%d',groundfl);
printFigFmt(h2,fig_path,fname,'eps');
%%%

% figure('Name','Field points');
% hold on
% imagesc(xMin:dx:xMax,zMin:dz:zMax,20*log10(abs(p_i_farXZ/pref)),...
%     'Interpolation','bilinear')
% contour(xMin:dx:xMax,zMin:dz:zMax,20*log10(abs(p_i_farXZ/pref)),...
%     10,'linecolor',[0 0 0])
% hold off
% xlabel('X')
% ylabel('Z')
% title('Reconstruction of the incident pressure over a plane grid XZ')
% axis equal
% colorbar
% colormap('jet')
% set(gca,'clim',[40 110])
