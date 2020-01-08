function [P_tot] = generateSimu(Rm,Rs,pp_simu)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Function simulating the pressure measured by an array of microphones. 
% Locations of microphones can be loaded from a text file (original
% microphone array) or calculated using a Matlab code (Not implemented
% yet). The calculation takes into account that the sphere is scattering
% the incident sound field.
%
% Input parameters :
%   - Rm : microphone locations arranged rowwise as [x1 y1 z1]
%   - Rs : source location (monopole)
%   - pp_simu : structure containing different parameters
%                
%       pp_simu = struct('MaxOrder',10,...
%                        'freq',100,...
%                        'SphereRadius',0.8,...
%                        'c',340);
% 
% see also scatteredPressure
%
% Théo Mariotte - 2019/10/21 - S-NAH (ENSIM)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
num_source = size(Rs,1);
num_mic = size(Rm,1);

P_crnt = zeros(num_mic,1);
P_tot = P_crnt;

x = Rm(:,1); y = Rm(:,2); z = Rm(:,3);
xs = Rs(:,1); ys = Rs(:,2); zs = Rs(:,3);

[rho_mic,phi_mic,theta_mic] = sphericalCoordinates(x,y,z);
[rho_src,phi_src,theta_src] = sphericalCoordinates(xs,ys,zs);

for isrc = 1 : num_source
    
    if pp_simu.incidentOnly
        
        r0 = [xs(isrc) ys(isrc) zs(isrc)];
        
        for imic = 1 : num_mic
            r = [x(imic) y(imic) z(imic)];
            P_crnt(imic) = incidentPressure(r,r0,pp_simu);
        end
        
    else

        % source location ins spherical coordinates
        r0 = [rho_src(isrc) theta_src(isrc) phi_src(isrc)];

        for imic = 1 : num_mic

            % microphone locations in spherical coordinates
            r = [rho_mic(imic) theta_mic(imic) phi_mic(imic)];

            % Sound pressure at the current microphone
            P_crnt(imic) = scatteredPressure(r,r0,pp_simu);       

        end

    end
    % Total sound field is the sum of each contibution (in the case of a
    % multisources simulation)
    
    P_tot = P_tot + P_crnt;
end


end