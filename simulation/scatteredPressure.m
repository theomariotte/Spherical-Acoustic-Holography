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
        
%         pp = struct(...
%             'pltype','3Dsphere',...
%             'valtype','real',...
%             'orders',[n m],...
%             'fontsize',12);
%         hh = SHvisualization(pp);
        
    end

    % spherical hankel function at the source location
    [hn_r0,~,~] = SphericalHankel1(n,k*r0_src);
    
    % derivative of the spherical hankel function on the surface of the
    % sphere
    [~,dhn_a,~] = SphericalHankel1(n,k*a);
    
    % first and second kind bessel function at the microphone location
    [jn_r,~,yn_r,~,~] = SphericalBessel(n,k*r_mic);
    
    % derivatives of Bessel functions on the surface of the sphere
    [~,djn_a,~,dyn_a,~] = SphericalBessel(n,k*a);
    
    
    switch pp_simu.method
        case 'GN'
             p_tmp = p_tmp + (hn_r0/dhn_a) * SH_sum;
             gain = (-1i*rho*c*Q)/(a^2);
             
        case 'Propagator'
            % propagator (see. Williams - Vector intensity recnstruction)
            Gn = (k*a)^2 * (jn_r * dyn_a - djn_a * yn_r);

            % sound pressure at the degree n
            p_tmp = p_tmp + Gn * ( hn_r0/dhn_a) * SH_sum;
            
            % gain that multiply the sound pressure
            gain = -4*pi / ((k*a)^2);
            
        case 'Brute'
            [hn_r,~,~] = SphericalHankel1(n,k*r_mic);
            Gn = ((jn_r * dhn_a - djn_a * hn_r)./dhn_a) .* hn_r0;
            p_tmp = p_tmp + Gn * SH_sum;
            gain = 4*pi*(1i)^n;
        otherwise
            warning('Default case : Green Neumann')
            p_tmp = hn_r0/dhn_a * SH_sum;
            gain = (-1i*rho*c*Q)/(a^2);
    end
   
end

% sound pressure measured by the current microphone
p = gain * p_tmp;

end