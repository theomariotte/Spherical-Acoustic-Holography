function [reconstructed_SF] = SNAH(theta_mic,phi_mic,P_meas,pp_reconstruction)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% [reconstructed_SF] = SNAH(theta_mic,phi_mic,P_meas,pp_reconstruction)
%
% Function allowing to reconstruct the sound field on a sphere of radius R
% given the measurements on a spherical microphone array of radius a.
% Here, the Spherical Nearfield Acoustic Holography (SNAH) method is used.
% Spherical harmonics are used to define the angular shape of the sound
% field. Spherical Bessel and Hankel functions are used to "retropropagate"
% the sound field for r > a. 
%
% Inputs :
%   - P_meas : array contaning the measured pressure at a given frequency
%   (prbem is solved in the frequency domain i.e. for stationary sound
%   sources).
%   - theta_mic : elevation angle of each microphone
%   - phi_mic : azimuthal angle for each microphone
%   - pp_reconstruction : reconstruction parameters and infos structure.
%   Parameters are defined below :
%       + pp_reconstruction.freq = working frequency
%       + pp_reconstruction.R = radius of the sphere where the
%       reconstruction has to be made (should be greater than a)
%       + pp_reconstruction.a = radius of the microphone array
%       + pp_reconstruction.c = sound wave celerity in the medium
%       + pp_reconstruction.maxOrder = maximum order N used to compute
%       spherical harmonics and Bessel functions.
%       + pp_reconstruction.incidentOnly = if 1, reconctruct only the
%       incident sound field (as if the array was transparent) ; if 0, the
%       total sound field is reconstructed.
%   
% Outputs
%   - P_reconstruct : reconstructed sound field on the sphere of radius
%   r = R. 
%
% see also sphericalTFSolver.m getSphericalHarmonics.m 
%
% Théo Mariotte - 11/2019 - ENSIM (SNAH project)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if nargin < 2
    error('Reconstruction parameters should be specified !');
end 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Define some parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Number of microphone on the sphere
Nmic = size(theta_mic,1);
% reconstruction radius
r_reconstruct = pp_reconstruction.R;
% radius of the sphere
a = pp_reconstruction.a;
% wave number
k = 2*pi*pp_reconstruction.freq / pp_reconstruction.c;
% maximal functions order
Nmax = pp_reconstruction.maxOrder;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Fourier Coefficients solving
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

idx = 1;
Ymat = zeros(Nmic,(Nmax+1)^2);
radial_fact = zeros((Nmax+1)^2,1);
idx_sup = 1;
for n = 0 : Nmax    
    % spherical hankel at the reconstruction radius
    [hn_r,~,~] = SphericalHankel2(n,k*r_reconstruct);        
    % derivative of the spherical hankel function on the surface of the
    % sphere
    [hn_a,dhn_a,~] = SphericalHankel2(n,k*a);    
    % first and second kind bessel function at the microphone location
    [jn_r,~,yn_r,~,~] = SphericalBessel(n,k*r_reconstruct);    
    % derivatives of Bessel functions on the surface of the sphere
    [jn_a,djn_a,~,~,~] = SphericalBessel(n,k*a);
    
%     fact_n = jn_r / (jn_a - (djn_a/dhn_a)*hn_a);
    if ~pp_reconstruction.incidentOnly
        fact_n = jn_r * dhn_a - djn_a*hn_r;
    else
        fact_n = jn_r * dhn_a;
    end
    
    idx_inf = idx_sup;
    idx_sup = idx_inf + 2*n;
    radial_fact(idx_inf : idx_sup) = repmat(fact_n,1,2*n+1);
    
   for m = -n : n
       % Matrice des harmoniques sphériques
        Y_tmp = getSphericalHarmonics(theta_mic,phi_mic,n,m); 
        Ymat(:,idx) = Y_tmp;
        
        % update
        idx = idx + 1;        
   end
   idx_sup = idx_sup+1;
end

% Compute Fourier coefficients for Spherical Fourier Transform (SFT)
Pmn = sphericalTFSolver(Ymat,P_meas);

% reconstruction
reconstructed_SF = 1j * (k*a)^2 * (Ymat * (radial_fact .* Pmn));

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%% Sound field reconstruction on a bigger sphere
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% % loop initialization
% idx_coef = 1;
% anglularSum = zeros(Nmic,1);
% P_tmp = anglularSum;
% 
% for n = 0 : Nmax
%     
%     % spherical hankel function at the source location
%     [hn_r,~,~] = SphericalHankel1(n,k*r_reconstruct);    
%     % derivative of the spherical hankel function on the surface of the
%     % sphere
%     [~,dhn_a,~] = SphericalHankel1(n,k*a);    
%     % first and second kind bessel function at the microphone location
%     [jn_r,~,~,~,~] = SphericalBessel(n,k*r_reconstruct);    
%     % derivatives of Bessel functions on the surface of the sphere
%     [jn_a,djn_a,~,~,~] = SphericalBessel(n,k*a);
%     
%     % propagator
%     if pp_reconstruction.incidentOnly
%         Gn = jn_r * dhn_a;
%     else
%         Gn = jn_r * dhn_a - djn_a * hn_r;
%     end
%     
%     for m = -n:n
%         anglularSum = anglularSum + Ymat(:,idx_coef) * Pmn(idx_coef);        
%         idx_coef = idx_coef + 1;
%     end
%     
%     P_tmp = P_tmp + Gn * anglularSum;
%     anglularSum = zeros(Nmic,1);
% end
% 
% reconstructed_SF = -1j * (k*a)^2 * P_tmp;

end