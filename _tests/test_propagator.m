% test propagator

clear; clc; 
close all

% paths
matlabhome = 'C:\Users\Théo\Documents\1_WORK\01_ENSIM\5A\Projet 5A\Spherical-Aoustic-Holography\';
simu_path = [matlabhome 'simulation\'];
plot_path = [matlabhome 'plot\'];
sph_func_path = [matlabhome 'spherical_functions\'];
data_path = [matlabhome 'data\'];

addpath(plot_path);
addpath(sph_func_path);
addpath(simu_path);
addpath(data_path);

a = 15e-2;

%% Test propagateur
r = 2*a;
k = linspace(0.001,100,2048);

n_vec = [4 7 9];

h = figure('Name','Propagator');
hold on
for i = 1 : length(n_vec)
    
    n = n_vec(i);
    
    % derivative of the spherical hankel function on the surface of the
    % sphere
    [~,dhn_a,~] = SphericalHankel1(n,k*a);

    % first and second kind bessel function at the microphone location
    [jn_r,~,yn_r,~,~] = SphericalBessel(n,k*r);

    % derivatives of Bessel functions on the surface of the sphere
    [~,djn_a,~,dyn_a,~] = SphericalBessel(n,k*a);

    % propagator
    Gn = (k*a).^2 .* (jn_r .* dyn_a - djn_a .* yn_r);
    
    logGn = 20*log10(abs(Gn));
    plot(k*a,logGn,'linewidth',2);    
    L{i} = sprintf('n = %d',n);
end
set(gca,'ylim',[-20 80])
grid on
xlabel('ka')
legend(L);

%% Test coefficients de Fourier analytique

theta0 = 0;
phi0 = 0;
r0 = 5*a;

ka = 2;

k = ka/a;
Nmax = 5;
n_vec = 0 : Nmax;

if theta0 == 0 && phi0 == 0
    Pmn = zeros(Nmax,1);
else
    Pmn = zeros(2*Nmax+1,Nmax);
end

u_vec = [0 1];

for ij = 1 : length(u_vec)
    
    u = u_vec(ij);
    fprintf('u = %d\n',u);
    
    for kk = 1 : length(n_vec)    

        % ordre courant
        n = n_vec(kk);
        m_vec = -n : n;

        % Hankel sphérique au point source
        [hn_kr0,~,~] = SphericalHankel1(n,k*r0);

        % Dérivée de Hankel sphérique sur la surface de la sphère
        [~,dhn_ka,~] = SphericalHankel1(n,k*a);

%         if theta0 == 0 && phi0 == 0
        if u == 0
            % cas de localisation sur l'axe polaire : tracé des coefficints en
            % focntion de la valeur de n
            Pmn(n+1) = -sqrt( 4*pi * (2*n+1) )*( hn_kr0 / (ka^2 * dhn_ka) );        
        else
            for m = 1 : length(m_vec)

                % harmonique sphérique au point source
                Ymn = getSphericalHarmonics(theta0,phi0,n,m_vec(m));

                % Coefficient de Fourier pour l'ordre n et le degré m (tous les
                % coefs d'ordre m ~= 0 doivent être nuls lorsque la source est sur
                % l'axe polaire (i.e. theta = phi = 0)
                Pmn(m,n+1) = -4*pi*(hn_kr0/(ka^2 * dhn_ka)) * Ymn;       

            end  

        end

    end
    
%     if theta0 == 0 && phi0 == 0
        figure
%         hold on
        if u == 0
            plot(n_vec,20*log10(abs(Pmn)),'k.','markersize',24)
        elseif u == 1
            Pmn_plt = diag(Pmn);
            plot(n_vec,20*log10(abs(Pmn_plt)),'k.','markersize',24)
        end
%         n_vec_int = 0:0.2:Nmax;
%         Pmn_int = interp1(n_vec,abs(Pmn),n_vec_int,'pchip');
%         plot(n_vec_int,20*log10(abs(Pmn_int)),'k','linewidth',2)
        grid on
        xlabel('n')
        ylabel('[dB]')
        set(gca,'ylim',[-60 0]);
%     end

end

