function [P_reconstruct] = SNAH(theta_m,phi_m,P_meas,pp_reconstruction)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function [P_reconstruct] = SNAH(P_meas,pp_reconstruction)
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
%   - theta_m : elevation angle of each microphone
%   - phi_m : azimuthal angle for each microphone
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
% see also solveIllConditionned.m getSphericalHarmonics.m 
%
% Théo Mariotte - 11/2019 - ENSIM (SNAH project)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if nargin < 2
    error('Reconstruction parameters should be specified !');
end 

%%% Define some parameters

% Number of microphone on the sphere
M = size(theta,1);
% reconstruction radius
R = pp_reconstruction.R;
% radius of the sphere
a = pp_reconstruction.a;
% wave number
k = 2*pi*pp.reconstruction.freq / pp_reconstruction.c;

% parameters for the Fourier coefficients computation. A pseudo inversion
% is used in the function
pp_fourier = struct('regularization',0,...
            'doplot',0,...
            'compute_condition',0);

% Get the spherical harmonics matrix Ymat (WIP)
for n = 1 : pp_reconstruction.maxOrder
   for m =  -n : n
       Ynm = getSphericalHarmonics(theta,phi,n,m);
       % Ymat has to be organized as follow
       %
       % Ymat = [[Y00(t1,p1) Y-11(t1,p1) ... YNN(t1,p1)] ;... 
       %         [Y00(t2,p2) Y-11(t2,p2) ... YNN(t2,p2) ;...
       %             :         :              :
       %             :         :              :
       %         [Y00(tM,pM) Y-11(tM,pM) ... YNN(tM,pM)]
   end
end

% Compute Fourier coefficients
Pmn = solveIllPosedProblem(Ymat,P_meas,pp_fourier);

% Compute the pressure on the new sphere using a propagator defined as 
%               Gn(a,r) = jn(kr)hn'(ka)-jn'(ka)hn(kr)
% See Jacobsen et al. 

% init
% Fourier coef index
idx_ = 1;
% reconstructed pressure on the sphere of radius R
P_reconstruct = zeros(size(theta));

for n = 1 : pp_reconstruction.maxOrder
    % Propagator 
    [jn_r,~,~,~,~] = SphericalBessel(n,k*R);
    [hn_r,~,~] = SphericalHankel1(n,k*R);
    [~,djn_a,~,~,~] = SphericalBessel(n,k*a);
    [~,dhn_a,~] = SphericalHankel1(n,k*a);
    Gn = jn_r*dhn_a - djn_a*hn_r;
        
    % init the sum over n
    S_tmp = zeros(size(theta));    
    for m =  -n : n
        % spherical harmonics at the current n & m for each microphone
        % location
        Ynm = getSphericalHarmonics(theta,phi,n,m);
        
        % update the sum
        S_tmp = S_tmp .*( Ymn * Pmn(idx_) );
        
        % upate Fourier coef index
        idx_ = idx_+1;
        
    end
    
    P_reconstruct = P_reconstruct + ( Gn * S_tmp );
    
end

P_reconstruct = -1i * (k*a)^2 * P_reconstruct;


end