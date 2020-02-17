function p = scatteredPressure(r,r0,pp_simu)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Function simulating the pressure measured by microphones flush mounted on
% a rigid sphere. 
% This measurement is simulated by computing the incident field and the
% pressure field scattered by the sphere. 
% Expression given by Williams in (DOI:10.1121/1.3278591).
%
% Théo Mariotte - 2019/10/21 - S-NAH (ENSIM)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% Initialization
p_tmp = 0;

%%% some parameters

% microphone array radius
a = pp_simu.SphereRadius;
% where the pressure is reconstructed
r_reconstruct = pp_simu.ReconstructRadius;

% wave number
f = pp_simu.freq;
w = 2*pi*f;
c = pp_simu.c;
rho = pp_simu.rho;
Q = pp_simu.Q;
k = w/c;

% microphone : 
% radial location
r_mic = r(1);
% Azimuth angle
theta_mic = r(2); 
% elevation angle
phi_mic = r(3); 

if (r_mic - r_reconstruct) > 1e-4
    error('Radial microphone location should be the same as the array radius !')
end

% monopole : 
% radial location
r0_src = r0(1);
% Azimuth angle
theta_src = r0(2); 
% elevation angle
phi_src = r0(3); 

for n = 0 : pp_simu.MaxOrder
    % sum of the product between spherical harmonics at the receiver
    % location and at the source location. 
    % The sum is reset to 0 at each loop over n (see Green Neumann)
    SH_sum = 0;
    
    for m = -n : n
       
        % compute the spherical harmonics for the current source and
        % receiver locations (theta = elevation ; phi = azimuth)
        Ynm_mic = getSphericalHarmonics(theta_mic,phi_mic,n,m);
        Ynm_src = getSphericalHarmonics(theta_src,phi_src,n,m);
        
        % compute the conjugate product and sum over all orders
        SH_sum = SH_sum + Ynm_mic * conj(Ynm_src);       
        
    end

    % spherical hankel function at the source location
    [hn_r0,~,~] = SphericalHankel2(n,k*r0_src);
    
    % derivative of the spherical hankel function on the surface of the
    % sphere
    [~,dhn_a,~] = SphericalHankel2(n,k*a);

    p_tmp = p_tmp + (hn_r0/dhn_a) * SH_sum;
             
   
end

% sound pressure measured by the current microphone
gain = -(1i*rho*w*Q)/(k*a^2);
p = gain * p_tmp;

end