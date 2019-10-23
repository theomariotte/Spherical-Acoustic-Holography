function [P] = genSimuMeas(xgrid,ygrid,zgrid,Rs,freq,SR_dist,pp_source)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Function simulating the sound pressure measured by a spherical array of
% microphones. 
% Input parameters : 
%   - Rm : microphones locations organized as [[x1 y1 z1];...;[xN yN zN]]
%   - Rs : sources locations (same organization as mic's locations)
%   - pp_source : physical qualities of each source (monopole)
%
% Théo Mariotte - 2019/10/18 - S-NAH (ENSIM)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Ns = size(Rs,1);
iMax = size(xgrid,1);
jMax = size(xgrid,2);
P = zeros(size(xgrid));
P_tmp = P;
Nmax = 10;

for s = 1 : Ns
    [phi_s,theta_s,r_s] = cart2sph(Rs(s,1),Rs(s,2),Rs(s,3));
    Rs_sph = [r_s,theta_s,phi_s];
    for ii = 1 : iMax
        for jj = 1 : jMax
                        
            switch pp_source.GreenFunc
                case 'GN'
                    [phi_m,theta_m,r_m] = cart2sph(xgrid(ii,jj),ygrid(ii,jj),zgrid(ii,jj));                      
                    Rm = [r_m,theta_m,phi_m];
                    P_tmp(ii,jj) = GreenNeumannFunction(Rs_sph,Rm,freq,Nmax,pp_source);
                case 'G'
                    r = [xgrid(ii,jj) ygrid(ii,jj) zgrid(ii,jj)];
                    P_tmp(ii,jj) = GreenFunction(r,Rs(s,:),freq,pp_source);
                otherwise
                    warning('Free field Green function used');
                    r = [xgrid(ii,jj) ygrid(ii,jj) zgrid(ii,jj)];
                    P_tmp(ii,jj) = GreenFunction(r,Rs(s,:),freq,pp_source);                    
            end     
            
        end
    end
    P = P + P_tmp;
end

end