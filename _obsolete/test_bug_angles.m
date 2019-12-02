% Test simulation mesure 
clear; clc; 
close all

% parameters
Nmax = 6;
f = 250;
Q = 1e-3;
c = 340;
k = 2*pi*f/c;

% source location
xs = 0; 
ys = 0;
zs = 1;

reconstruct = 1;
r_reconstruct = 30e-2;
typ = 1;

%% initialisation

if typ == 1
    % ouverture fichier
    Nmic = 36;
    data_path = 'data/';
    filename =  '3Dcam36Officiel.txt';
    fID = fopen([data_path filename],'r');
    sizeA = [Nmic 3];
    mic_loc_tmp = fscanf(fID,'%f');
    fclose(fID);
    
    Rm = transpose(reshape(mic_loc_tmp,[3 Nmic]));  
    dointerp = 1;
    
    figure
    plot3(Rm(:,1),Rm(:,2),Rm(:,3),'kx')
    xlabel('X')
    ylabel('Y')
    zlabel('Z')
    grid on
    axis equal

else
    
    ntheta = 12;
    nphi = 11;
    doplot = 1;
    a = 15e-2;

    Rm = buildSimuSphere(ntheta,nphi,a,doplot);
    
    Nmic = ntheta * nphi;
    
end

x = Rm(:,1);
y = Rm(:,2);
z = Rm(:,3);

[rho,phi,theta] = sphericalCoordinates(x,y,z);
[rhos,phis,thetas] = sphericalCoordinates(xs,ys,zs);

% structure
a = rho(1);
pp_simu = struct('MaxOrder',Nmax,...
                'freq',f,...
                'SphereRadius',a,...
                'c',c,...
                'rho',1.2,...
                'Q',Q,...
                'method','Propagator',...
                'doplot',0);

%% simulation
r = [rho theta phi];
r0 = [rhos thetas phis];
for imic = 1 : Nmic
    r_crnt = r(imic,:);
    P_crnt(imic) = scatteredPressure(r_crnt,r0,pp_simu);
end
  
param = struct('sphereRadius',a,...
            'numAngle',30,...
            'interpType','natural',...
            'doplot',1,...
            'theta0',thetas,...
            'phi0',phis);

[P_interp,x_grid,y_grid,z_grid] = ...
    fieldInterpolation(real(P_crnt'),Rm,param);

%% Coefficients de Fourier
idx = 1;
Ymat = zeros(Nmic,(Nmax+1)^2);
Pmn_th = zeros((Nmax+1)^2,1);

for n = 0 : Nmax
    
    % pour Fourier
    [hn_r0,~,~] = SphericalHankel1(n,k*rhos);
    [~,dhn_a,~] = SphericalHankel1(n,k*a);
    
   for m = -n : n
       % Matrice des harmoniques sphériques
        Y_tmp = getSphericalHarmonics(theta,phi,n,m); 
        Ymat(:,idx) = Y_tmp;
                
        % coefficients de Fourier analytiques
        Y_source = getSphericalHarmonics(thetas,phis,n,m);
        Pmn_th(idx) = -4*pi*hn_r0/((k*a)^2*dhn_a) * conj(Y_source); 
        
        % update
        idx = idx + 1;        
   end
end
% 
% fourier coefficients
pp = struct('regularization',0,...
           'doplot',0,...
           'compute_condition',1);
[Pmn,cond_num] = solveIllPosedProblem(Ymat ,P_crnt,pp);
% 
% %% reconstruction
% idx = 1;
% P_reconstruct = 0;
% 
% for n = 0 : Nmax
%     
%     [jn_r,~,yn_r,~,~] = SphericalBessel(n,k*r_reconstruct);
%     [~,djn_a,~,dyn_a,~] = SphericalBessel(n,k*a);
%     [hn_r,~,~] = SphericalHankel1(n,k*r_reconstruct);
%     [~,dhn_a,~] = SphericalHankel1(n,k*a);
%     
%     if reconstruct == 1
%         Gn = -1j*(k*a)^2 * (jn_r * dhn_a - djn_a*hn_r);
%     else
%         Gn = (k*a)^2 * (jn_r * dyn_a - djn_a*yn_r);
%     end
%     
%     S_tmp = 0;
%     for m = -n : n                
%         S_tmp = S_tmp + ( Ymat(:,idx) * Pmn_th(idx) );
%         idx = idx + 1;
%     end
%     
%     P_reconstruct = P_reconstruct + S_tmp;
%     
% end
% 
% % interpolation sur une nouvelle grille pour l'affichage
% [rho_tmp,phi_tmp,theta_tmp] = cart2sph(Rm(:,1),Rm(:,2),Rm(:,3));
% 
% [x_r,y_r,z_r] = sph2cart(r_reconstruct,phi_tmp,theta_tmp);
% param = struct('sphereRadius',r_reconstruct,...
%             'numAngle',70,...
%             'interpType','natural');            
% [P_r_interp,x_grid_r,y_grid_r,z_grid_r] = fieldInterpolation(P_reconstruct,[x_r y_r z_r],param);

%% figures

figure('Name','Sound pressure on the sphere');
hold on
surf(x_grid,y_grid,z_grid,20*log10(abs(P_interp)/20e-6))    
% surf(x_grid,y_grid,z_grid,real(P_interp))    
colormap('jet')
shading('interp')
cc = colorbar;
cc.Label.String = 'Sound Pressure Level [dB]';
axis equal
title('Measured sound pressure')
plot3(xs,ys,zs,'ro','linewidth',3)
delta = 0.1;
plot3(x+delta,y+delta,z+delta,'k.','linewidth',5,'markersize',12);
hold off
xlabel('X','interpreter','latex')
ylabel('Y','interpreter','latex')
zlabel('Z','interpreter','latex')
axis('equal')
grid on

% % tracé des coefficients
nn = linspace(0,Nmax,(Nmax+1)^2);
hh = figure('Name','Fourier coefficient');
hold on
    plot(20*log10(abs(Pmn_th(Pmn_th ~=0))),'kx');
    plot(20*log10(abs(Pmn(Pmn ~=0))),'ro');    
hold off
xlabel('n')
ylabel('20log_{10}(|P_{mn}|)')
title(sprintf('Fourier coefficients : ka = %2.1f',k*a));
legend('Theoretical','Computed')
grid on

% 
% % figures
% figure('Name','Reconstructed sound pressure');
% hold on
% surf(x_grid_r,y_grid_r,z_grid_r,20*log10(abs(P_r_interp)/20e-6))    
% colormap('jet')
% shading('interp')
% cc = colorbar;
% cc.Label.String = 'Sound Pressure Level [dB]';
% axis equal
% title(sprintf('Reconstructed sound pressure : r = %d',r_reconstruct))
% 
% plot3(xs,ys,zs,'ro','linewidth',3)
% hold off
% xlabel('X','interpreter','latex')
% ylabel('Y','interpreter','latex')
% zlabel('Z','interpreter','latex')
% axis('equal')
% grid on
% hold off

% figure
% plot(theta,real(P_crnt),'kx','linewidth',2)
% xlabel('Elevation angle $\theta$','interpreter','latex')
% ylabel('Sound pressure $|p(\mathbf(r))|$','interpreter','latex')
% 
% figure
% plot(phi,real(P_crnt),'kx','linewidth',2)
% xlabel('Azimuthal angle $\phi$','interpreter','latex')
% ylabel('Sound pressure $|p(\mathbf(r))|$','interpreter','latex')




