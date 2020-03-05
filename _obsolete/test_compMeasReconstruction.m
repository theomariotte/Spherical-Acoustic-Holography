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

%% Pressure measurements over reconstruction axis

% reference sound pressure
dataPath = ['/Users/theomariotte/Documents/01_work/ENSIM/5A/Projet/'...
    'Spherical-Acoustic-Holography/data/'];


no_devices=3; % n°carte
no_voies{1}=5; % n°voie
nbvoies2=1;

nameVec = {'12022020_112037',...
'12022020_112239',...
'12022020_112443',...
'12022020_112548',...
'12022020_112648',...
'12022020_112732',...
'12022020_112909',...
'12022020_113046',...
'12022020_113211',...
'12022020_113318'};

for ii = 1 : length(nameVec)
    cd([dataPath nameVec{ii} '/'])
    [pRefTime,tref,~]=convert_data_TDMS(nbvoies2,type_data,no_devices,no_voies,0);
    pMeasAxis(:,ii) = shiftdim(pRefTime.voie1);
    pause(0.001);
end
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

for ii = 1 : nbvoies
   
    crnt_x = eval(sprintf('pArrayTime.voie%d',ii));
    Pxy(:,ii) = cpsd(pRefTime.voie1,crnt_x,hann(M),M/2,M);
    pMeas(ii) = Pxy(idxF,ii);
    %plot(ff,10*log10(abs(Pxy(:,k)/(pref^2))));
end

% grid on
% set(gca,'xlim',[0 fs/2]);

%% Spectrum over reconstruction axis
nbFieldPoints = size(pMeasAxis,2);

Pxy = zeros(M/2+1,nbFieldPoints);
pMeasAxisFreq = zeros(nbFieldPoints,1);

ff = linspace(0,fs/2,M/2+1);
[~,idxF] = min(abs(ff-f0));

% figure('Name','PSD')
% hold on

for ii = 1 : nbFieldPoints
   
    crnt_x = eval(sprintf('pArrayTime.voie%d',ii));
    Pxy(:,ii) = cpsd(pRefTime.voie1,crnt_x,hann(M),M/2,M);
    pMeasAxisFreq(ii) = Pxy(idxF,ii);
    %plot(ff,10*log10(abs(Pxy(:,k)/(pref^2))));
end


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

%%% transfert function between sources and array 

% green neumann without taking the ground into account
GN = getGreenNeumannFRF(ES_loc,mic_loc,Nmax,k,a);
    
% reg param = 1e-2 good with tikhonov !
q = solveIllPosedProblem(GN,pMeas,'tikhonov',5e-3);


%% reconstruction over the reconstruction axis

% axis where the point measurement was done
dMin = 5e-2;
dMax = 50e-2;
dy = 5e-2;

expAxis = ys - (dMin : dy : dMax);

% axis where the sound field has to be reconstructed
dy = 1e-3;
y_reconstruc = (ys - (dMin : dy : dMax))';

x_reconstruc = zeros(size(y_reconstruc));
z_reconstruc = zeros(size(y_reconstruc));

R_reconstruc = [x_reconstruc y_reconstruc z_reconstruc];

G = getFreeFieldFRF(ES_loc_cart,R_reconstruc,k);

% incident sound field reconstruction
p_i = G * q;


%% Figures

pp_interp = struct('sphereRadius',a,...
            'numAngle',30,...
            'interpType','natural',...
            'doplot',0);

pp_plot = struct('showSource',0,...
                'showMic',0,...
                'plt_typ','dB',...
                'clim',[40 65],...
                'fontSize',12,...
                'freq',f0);              

handle1 = pressureMeasurementVisu(pMeas,Rm,Rs,pp_interp,pp_plot);
colormap('jet')
view(135,17);
xlabel('Y')
ylabel('X')
zlabel('Z')
fname = sprintf('measuredSoundField%dHz',f0);
printFigFmt(handle1,fig_path,fname,'png')

% hold(gca,'on')
% plot3(X_es,Y_es,Z_es,'o','color',[.5 .5 .5]);
% if groundfl
%    plot3(IES_loc_cart(:,1),IES_loc_cart(:,2),IES_loc_cart(:,3),...
%        'o','color',[.5 .5 .5]); 
% end
% [xGround,yGround] = meshgrid(-1:.5:1);
% zGround = zeros(size(xGround)) - H;
% surf(xGround,yGround,zGround);
%hold(gca,'off')

%%%%%


q_mat = reshape(q,sz);

handle2 = figure('Name','Equivalent sources strength');
hold on
imagesc(Y_es_tmp,Z_es_tmp,abs(q_mat),'Interpolation','bilinear')
set(gca,'yDir','normal')
dessin_13HP(0,0.05)
hold off
colormap jet
xlabel('Y')
ylabel('Z')
axis equal
set(gca,'fontsize',12,'xlim',[min(Y_es_tmp) max(Y_es_tmp)],'ylim',[min(Z_es_tmp) max(Z_es_tmp)]);
fname = sprintf('EqSourcesStrength%dHz',f0);
printFigFmt(handle2,fig_path,fname,'eps');
%%%%

% hh = figure('Name','Sound pressure over reconstruction axis');
% hold on
% plot(expAxis,20*log10(abs(pMeasAxisFreq/pref)),'--k','linewidth',2,'marker','o')
% plot(y_reconstruc,20*log10(abs(p_i/pref)),'color',[.5 .5 .5],'linewidth',3);
% hold off
% xlabel('Distance to the source [m]','interpreter','latex')
% ylabel('SPL [dB]','interpreter','latex')
% grid on
% legend('Simulation','Reconstruction')
% title('Reconstruction using LSQR')
set(gca,'fontsize',12)
% fname = sprintf('reconstructionXaxisRSB%d',round(SNR));
% printFigFmt(hh,fig_path,fname,'.eps');
