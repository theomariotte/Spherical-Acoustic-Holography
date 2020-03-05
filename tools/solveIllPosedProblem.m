function [q] = solveIllPosedProblem(G,p,regType,regParam)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% [Pmn] = solveIllPosedProblem(X,s,regType,regParam)
%
% This function solves an ill-posed problem. That means that the system 
%                               G.q = p
% cannot be solved as a linear system. 
%
% To solve this problem, three regularization methods are investigated : 
%   - Tickhonov which consists in solving G' (G.G' + uI)^-1 p, where "u" is
%   the regularization parameter 
%
%   - lsqr which sole Gq = p using iterative algorithm (based on conjugate
%   gradient)
%
%   - TSVD which truncates the singular value decomposition to r values.
%   Here, r is the regularization parameter.
%
% Input :
%   * G : FRF matrix to be inverted
%   * p : measured pressure vector
%   * regType : type of regularization ('tikhonov','lsqr','TSVD')
%   * regParam : regularization parameter (should match with the choosen
%   regularization method)
%
% Output :
%   * q : solution of the problem
%
% Théo Mariotte - 2019/10/24 - S-NAH (ENSIM)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% Initialization
if nargin < 4
    % If no regularization parameter is given, simple pseudo-inversion is
    % computed
    if nargin == 3
        switch regType
            case 'tikhonov'
                regParam = 1e-5;
            case 'lsqr'
                regParam = 1e-6;
            case 'tsvd'
                regParam = 10;
            otherwise
                regType = 'lsqr';
                regParam = 1e-6;
                warning('Default LSQR method used with tol = %1.2e',regParam);
        end
    else
        error('Please define a regularization method');
    end
    
    warning('Reg method %s : default reg. param used r=%1.2g',regType,regParam);
end

if nargin < 2
    error('One matrix and one vector should be defined to solve the problem !')
end

p = shiftdim(p);

if size(p,2) > 1,error('Data vector p must be 1-D array');end

if size(G,1) ~= size(p,1)
    if size(G,2) == size(p,1)
        G = G';
    else
       error('Number of rows in G should be the same as in p') 
    end
end
 
switch regType
    
    case 'tikhonov'
        
        G_trans = G';
        G_prod = G * G_trans;        

        I = eye(size(G_prod));

        G_whole = (G_prod + regParam * I);

        % inverse
        X_inv = pinv(G_whole) ;
        q = G' * (X_inv * p);
        
    case 'lsqr'
        
        nmax = 200;
        q = lsqr(G,p,regParam,nmax);
        
    case 'tsvd'
        
        [U,S,V] = svd(G);
        
        if regParam < 1, error('TSVD : regularization parameter should be a positive integer !');end
        
        % if the reg parameter is not an integer value
        r = round(regParam);
        
        % truncation of sigular values 
        U = U(:,1:r);
        V = V(:,1:r);
        S = S(1:r,1:r);
        
        X_inv = V * ( 1./S * U' );
        
        q = X_inv * p;
        
    otherwise
        regParam = 1e-6;
        nmax = 200;
        q = lsqr(G,p,regParam,nmax);

        warning('No regularization method found. LSQR as default with %1.2e tolerance.',regParam);
        
end

end