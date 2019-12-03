function [q,GN] = getSourceStrength(P_meas,mic_loc,src_loc,pp)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% [q] = getSourceStrength(P_meas,mic_loc,src_loc,pp)
%
% Function allowing to compute the 'strength' of equivalent source when an
% spherical equivalent sources model is solved for sound feld
% reconstruction. 
%
% input : 
%   - P_meas : column arranged vector containing measured sound pressure at
%   each microphone location on the sphere, in the frequency domain. 
%   - mic_loc : location of each microphone int spherical coordinates.
%   These data should be arranged linewise. Each column is associated to
%   one coordinate (rho, theta or phi). 
%   - src_loc : equivalent sources locations. Arranged in the same way as
%   mic_loc
%   - pp : structure with some parameters (see default)
%
% Output : 
%   - q : equivalent sources strength
%
% Théo Mariotte - 2019/12/03 - ENSIM (SNAH Project)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Initialization
M = size(mic_loc,1);
L = size(src_loc,1);

% build the Green Neumann function matrix
GN = zeros(M,L);
for l = 1 : L
    for s = 1 : M
        GN(s,l) = GreenNeumannFunction(src_loc(l,:),mic_loc(s,:),pp); 
    end
end

% ill posed system solver
[q,cdn] = solveIllPosedProblem(GN,P_meas,pp.reg_parameter);

fprintf('Ill-posed system condition number : %.2f\n',cdn)

end