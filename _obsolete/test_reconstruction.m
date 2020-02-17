% test simulation mesures
clear; clc;
close all;

%% Paramètres

data_path = 'data/';
fig_path = '';
% Tracé de la pression acoustique
% plt_typ : ('real') real part ; ('dB') decibels ; ('abs') absolute value
plt_typ = 'dB';

% Méthode de simulation
% ('GN') : utilisation de la fonction de Green Neumann (calcul possible
% uniquement sur la frontière) ; 
% ('Propagator') Utilise le propagateur
% défini par Williams dans 'Intensity vector reconstruction'.
% cete dernière permet de simuler la pression acoustique sur tout le
% domaine de valdité, contrairement à la précédente qui ne permet de
% calculer la pression que sur une sphère de rayon a
sim_method = 'GN';

% Microphone array radius
a = 15e-2;
ka = 2;

% Maximum order (Bessel, Hankel, SH)
Nmax = 4;
% Noise amp
nAmp = 0.00;
% color limits
clims = [50 90];
% reconstruction incident ou non
incidentOnly = 1;

%%% Autres paramètres
% débit de la source [m^3/s]
Q = 5e-4;
% masse volumique du milieu
rho_air = 1.2;
% célérité du milieu 
c = 340;
% fréquence de travail [Hz]
f = ka*c/(2*pi*a);

%% Chargement ds coordonnées des microphones (cartésien)

Nmic = 36;
filename =  '3Dcam36Officiel.txt';
fID = fopen([data_path filename],'r');
sizeA = [Nmic 3];
mic_loc_tmp = fscanf(fID,'%f');
fclose(fID);

Rm = transpose(reshape(mic_loc_tmp,[3 Nmic]));  

[rho,phi,theta] = sphericalCoordinates(Rm(:,1),Rm(:,2),Rm(:,3));

%% Source location

front_of_mic = 0;

if front_of_mic
    rhos = 0.8;
    phis = phi(1);
    thetas = theta(1);

    [xs,ys,zs] = cartesianCoordinates(rhos,phis,thetas);
    Rs = [xs ys zs];
else
    % cartesian coordinates [m]
    xs = 1;
    ys = [0];
    zs = [0];

    Rs = [xs ys zs];
    [rhos,phis,thetas] = sphericalCoordinates(xs,ys,zs);
end


%% calculate the soundfield at each microphone
dx = 1e-2;
r_reconstruct_vec = a : dx : rhos-dx;
% r_reconstruct_vec = 2*a;

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

% RSB computation
num_RSB = max(max(abs(P)));
RSB = 20*log10(num_RSB/nAmp);

fprintf('RSB = %2.0f dB\n',RSB)

P_meas = P + nAmp * randn(size(P));

%% reconstruction du champ acoustique 

N_reconstruct = length(r_reconstruct_vec);
reconstruct_SF = zeros(Nmic,N_reconstruct);

for k = 1 : N_reconstruct
    r_reconstruct = r_reconstruct_vec(k);
    pp_reconstruction = struct('a',a,...
                            'R',r_reconstruct,...
                            'freq',f,...
                            'c',c,...
                            'maxOrder',Nmax,...
                            'incidentOnly',incidentOnly);
    P_tmp = SNAH(theta,phi,P_meas,pp_reconstruction);
    reconstruct_SF(:,k) = P_tmp;
end

%% Figures

% pression simulée sur la sphère (interpolée pour une lecture plus simple)
pp_interp = struct('sphereRadius',a,...
            'numAngle',30,...
            'interpType','natural',...
            'doplot',0);
        
pp_plot = struct('showSource',1,...
                'showMic',1,...
                'plt_typ',plt_typ,...
                'clim',clims,...
                'fontSize',12,...
                'freq',f);     
            
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

handle1 = pressureMeasurementVisu(P_meas,Rm,Rs,pp_interp,pp_plot);

%%%
% Cas où l'on ne reconstruit qu'à 1 seul rayon
%%%
if size(reconstruct_SF,2) == 1
    % Pression reconstruite par méthode d'holographie acoustique sphérique
    pp_interp = struct('sphereRadius',r_reconstruct,...
                'numAngle',30,...
                'interpType','natural',...
                'doplot',0);

    pp_plot = struct('showSource',1,...
                    'showMic',0,...
                    'plt_typ',plt_typ,...
                    'clim',clims,...
                    'fontSize',12,...
                    'freq',f);     

    handle2 = pressureMeasurementVisu(reconstruct_SF,Rm,Rs,pp_interp,pp_plot);
    title('Recontructed Sound Field');
    colormap('jet')
    view(-20,40);
    
    % Pression reconstruite attendue (calculée par simulation)
    pp_simu = struct('MaxOrder',Nmax,...
                    'freq',f,...
                    'SphereRadius',a,...
                    'ReconstructRadius',r_reconstruct,...
                    'c',c,...
                    'rho',rho_air,...
                    'Q',Q,...
                    'method',sim_method,...
                    'incidentOnly',incidentOnly,...
                    'doplot',0);
    % simulation de la pression reconstruite
    [P_cmp] = generateSimu(Rm,Rs,pp_simu);  
    
    % interpolation du champ
    pp_interp = struct('sphereRadius',r_reconstruct,...
                'numAngle',30,...
                'interpType','natural',...
                'doplot',0);
    
   % Affichage du champ
    pp_plot = struct('showSource',1,...
                    'showMic',0,...
                    'plt_typ',plt_typ,...
                    'clim',clims,...
                    'fontSize',12,...
                    'freq',f);     

    handle3 = pressureMeasurementVisu(P_cmp,Rm,Rs,pp_interp,pp_plot);
    title('Expected Recontructed Sound Field');
    colormap('jet')
    view(-20,40);
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%
% Cas où la reconstruction est calculée sur plusieurs sphères de rayon a.
% Cela permet de comparer la reconstruction avec la valeur du champ
% attendue sur un point donné de la sphère de microphones 
%%%
else
    
    % point face à la source
    [~,idx_phi] = min(abs(phi - phis)) ;
    [~,idx_theta] = min(abs(theta - thetas)) ;
    
    % pression théorique
%     h = figure('Name','Evolution reconstruction');
    for k = 1 : N_reconstruct
        r_reconstruct = r_reconstruct_vec(k);
        
        % pression incidente reconstruite et interpolée        
        [P_interp,x_int,y_int,z_int] = fieldInterpolation(reconstruct_SF(:,k),Rm,pp_interp);
        
        %%%
        % On ne reconstruit que le champ incident
        %%%
        if incidentOnly
            % Pression incidente théorique        
            [x_grid,y_grid,z_grid] = cartesianCoordinates(r_reconstruct,phi(idx_phi),theta(idx_phi));
            Rm_grid = [x_grid,y_grid,z_grid];
            R = norm(Rm_grid - Rs);
            w = (2*pi*f);
            wave_num = w/c;
            switch sim_method
                case 'GN'
                    P_reconstruct_expect = (1i*w*rho_air*Q) * exp(-1i*wave_num*R)/(4*pi*R);
                case 'Propagator'
                    P_reconstruct_expect = exp(-1i*wave_num*R)/(4*pi*R);

            end
            P_reconstruct(k) = max(max(abs(P_interp)));
            
        %%%
        % Le champ total est reconstuit (incident + diffusé)
        %%%            
        else
            %simulation pour comparer        
            pp_simu = struct('MaxOrder',Nmax,...
                        'freq',f,...
                        'SphereRadius',a,...
                        'ReconstructRadius',r_reconstruct,...
                        'c',c,...
                        'rho',rho_air,...
                        'Q',Q,...
                        'method',sim_method,...
                        'incidentOnly',incidentOnly,...
                        'doplot',0);
            [P_cmp] = generateSimu(Rm,Rs,pp_simu);
            [P_interp_exp,x_int,y_int,z_int] = fieldInterpolation(P_cmp,Rm,pp_interp);

            P_reconstruct(k) = max(max(abs(P_interp)));
            P_reconstruct_expect(k) = max(max(abs(P_interp_exp)));
        end

%        % Pression reconstruite par méthode d'holographie acoustique sphérique
%         pp_interp = struct('sphereRadius',r_reconstruct,...
%                     'numAngle',30,...
%                     'interpType','natural',...
%                     'doplot',0);
% 
%         pp_plot = struct('showSource',1,...
%                         'showMic',0,...
%                         'plt_typ',plt_typ,...
%                         'clim',clims,...
%                         'fontSize',12,...
%                         'freq',f); 
%                     
%         
%         pressureMeasurementVisu(reconstruct_SF(:,k),Rm,Rs,pp_interp,pp_plot);
%         
%         hold on
%         [xx,yy,zz] = cartesianCoordinates(r_reconstruct,phi(idx_phi),theta(idx_phi));
%         plot3(xx,yy,zz,'bo','markersize',12,'linewidth',2)
%         hold off
%         title(sprintf('Recontructed Sound Field r = %.1e',r_reconstruct));
%         pause(.1)       
        
    end
  
    
    
    hh = figure('Name','Sound field reconstruction along the source axis');
    pref = 20e-6;
%     P_dB = 20*log10(abs(reconstruct_SF(idx_phi,:))/pref);
    P_dB = 20*log10(abs(P_reconstruct)/pref);
%     Pi_dB = 20*log10(abs(P_i)/pref);
    Pi_dB = 20*log10(abs(P_reconstruct_expect)/pref);
    
    r_source = fliplr(rhos - r_reconstruct_vec);
    hold on
    plot(r_source,fliplr(P_dB),'color',[.5 .5 .5],'linewidth',3)
    plot(r_source,fliplr(Pi_dB),'--k','linewidth',2)
    hold off
    grid on
    xlabel('Distance to the source [m]')
    ylabel('SPL [dB]')
    title(sprintf('ka = %1.1f ; N = %d',ka,Nmax))
    set(gca,'fontsize',12)
    legend('Reconstruction','Simulation')
%     fig_path = 'C:\Users\Théo\Documents\1_WORK\01_ENSIM\5A\Projet 5A\Rapports\Fiches de suivi\1912_breve_avancement\';
%     fname = sprintf('SNAH_f_%d',round(f));
%     printFigFmt(hh,fig_path,fname,'eps');

end
