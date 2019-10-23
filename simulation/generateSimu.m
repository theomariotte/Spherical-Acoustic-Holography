function [P] = generateSimu(Rm,Rs,pp_simu)
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
P = P_crnt;

figure('name','verif points')
hold on

for isrc = 1 : num_source
            
    [phi_s,theta_s,r_s] = cart2sph(Rs(isrc,1),Rs(isrc,2),Rs(isrc,3));
    r0 = [r_s theta_s phi_s];
    
    if pp_simu.doplot
        [xs_plt,ys_plt,zs_plt] = sph2cart(phi_s,theta_s,r_s);    
        plot3(xs_plt,ys_plt,zs_plt,'gx')
        plot3(Rs(isrc,1),Rs(isrc,2),Rs(isrc,3),'ro','linewidth',3)
        axis equal
        set(gca,'xlim',[-0.2 0.2],'ylim',[-0.2 0.2],'zlim',[-0.2 0.2])
        xlabel('X')
        ylabel('Y')
        zlabel('Z')
        title('Vérification points de chaque mesure')
    end
    
    for imic = 1 : num_mic
        
        [phi_m,theta_m,r_m] = cart2sph(Rm(imic,1),Rm(imic,2),Rm(imic,3));
        r = [r_m theta_m phi_m];
        
        % Sound pressure at the current microphone
        P_crnt(imic) = scatteredPressure(r,r0,pp_simu);
        
        if pp_simu.doplot
            [xr_plt,yr_plt,zr_plt] = sph2cart(phi_m,theta_m,r_m);    
            plot3(xr_plt,yr_plt,zr_plt,'rx')
            plot3(Rm(imic,1),Rm(imic,2),Rm(imic,3),'ko')
            grid on
            view(30 + imic * 0.5,30);            
            pause(.1)
        end
        
    end
    
    P = P + P_crnt;
end


end