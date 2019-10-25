function [varargout] = SphericalFourierCoeff(Ymn,p_vec,pp)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% [Pmn] = SphericalFourierCoeff(Ymn,p_vec)
%
% Function calculating the (N+1)^2 Spherical Fourier coefficients Pmn from 
% a spherical harmonics matrix Ymn  [N x (N+1)^2] and a pressure array of N
% measured values. 
% The problem is ill-posed and required a pseudo-inversion based on
% singular values decomposition.
%
% see also pinv
%
% Théo Mariotte - 2019/10/24 - S-NAH (ENSIM)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% Initialization
if nargin < 3
   pp = struct('regularization',0,...
               'doplot',1,...
               'compute_condition',1);
end

if nargin < 2
    error('Acouctic pressure and Spherical Harmonics should be defined !')
end

varargout = zeros(nargout,1);

[M,tst] = size(p_vec);

if tst > M
    p_vec = p_vec';
    [M,tst] = size(p_vec);
end

if tst > 1
    error('Acoustic pressure must be 1-D array')
end

[M2,L] = size(Ymn);

if M2 > L
    Ymn = Ymn';
    [M2,L] = size(Ymn);
end

if M2 ~= M, error('Ymn : number of lines should be the same as p_vec'); end

%%% Compute

if pp.regularization 
    
    Y_prod = Y * (conj(Y))';
    I = eye(size(Y_prod));
    
    Ymn_reg = Y_prod + alpha * I;
    
    if pp.compute_condition
        [~,S,~] = svd(Ymn_reg);
        max_S = max(max(S));
        min_S = min(min(S));
        condition = max_S/min_S;
        varargout(2) = condition;
    end
    
    Ymn_reg_inv = pinv(Ymn_reg);
    
    Ymn_inv = (conj(Y))' * Ymn_reg_inv ;
    
else
    
    if pp.compute_condition
        [~,S,~] = svd(Ymn);
        max_S = max(max(S));
        min_S = min(min(S));
        condition = max_S/min_S;
        varargout(2) = condition;
    end
    
    Ymn_inv = pinv(Ymn);    
        
end

Pmn = Ymn_inv * p_vec;
varargout(1) = Pmn;

%%% plot

if pp.doplot
   h = figure('Name','Spherical Fourier Coefficients');
   plot(Pmn,'k.','markersize',24)
   grid on
   title('Fourier coefficients $P_{mn}','interpreter','latex')   
end

end