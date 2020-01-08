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
SH_sum = 0;
p_tmp = 0;

%%% some parameters

% microphone array radius
a = pp_simu.SphereRadius;
% where the pressure is reconstructed
r_reconstruct = pp_simu.ReconstructRadius;

% wave number
f = pp_simu.freq;
c = pp_simu.c;
rho = pp_simu.rho;
Q = pp_simu.Q;
k = (2*pi*f)/c;

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
   
    for m = -n : n
       
        % compute the spherical harmonics for the current source and
        % receiver locations (theta = elevation ; phi = azimuth)
        Ynm_mic = getSphericalHarmonics(theta_mic,phi_mic,n,m);
        Ynm_src = getSphericalHarmonics(theta_src,phi_src,n,m);
        
        % compute the conjugate product and sum over all orders
        SH_sum = SH_sum + Ynm_mic * conj(Ynm_src);       
        
    end

    % spherical hankel function at the source location
    [hn_r0,~,~] = SphericalHankel1(n,k*r0_src);
    
    % spherical hnkel function at the receiver location
    [hn_r,~,~] = SphericalHankel1(n,k*r_reconstruct);
    
    % derivative of the spherical hankel function on the surface of the
    % sphere
    [~,dhn_a,~] = SphericalHankel1(n,k*a);
    
    % first and second kind bessel function at the microphone location
    [jn_r,~,yn_r,~,~] = SphericalBessel(n,k*r_reconstruct);
    
    % derivatives of Bessel functions on the surface of the sphere
    [~,djn_a,~,dyn_a,~] = SphericalBessel(n,k*a);
    
    
    switch pp_simu.method
        case 'GN'
             p_tmp = p_tmp + (conj(hn_r0)/conj(dhn_a)) * SH_sum;
             gain = -(1i*rho*c*Q)/(a^2);
             
        case 'Propagator'
            % cas ou r = a
            if r_reconstruct == a
                
                Gn = hn_r0/dhn_a;
                
                p_tmp = p_tmp + Gn * SH_sum;
                
                gain = -4*pi/(k*a)^2;
                
            % cas ou r != a
            else
                % propagator (see. Williams - Vector intensity recnstruction)
                Gn = (jn_r * dyn_a - djn_a * yn_r);
%                 Gn = (jn_r - (djn_a/dhn_a)*hn_r) * hn_r0;
                
                % sound pressure at the degree n
                p_tmp = p_tmp + Gn * ( hn_r0/dhn_a ) * SH_sum;
%                 p_tmp = p_tmp + Gn * SH_sum;

                % gain that multiply the sound pressure
                gain = -4*pi;
            end

            
        otherwise
            warning('Default case : Green Neumann')
             p_tmp = p_tmp + (conj(hn_r0)/conj(dhn_a)) * SH_sum;
             gain = -(1i*rho*c*Q)/(a^2);
    end
   
end

% sound pressure measured by the current microphone
p = gain * p_tmp;

end