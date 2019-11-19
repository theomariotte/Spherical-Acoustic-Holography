function [q,varargout] = solveIllPosedProblem(X,s,pp)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% [Pmn] = solveIllPosedProblem(Ymn,p_vec)
%
% This function solves an ill-posed problem. That means that the system 
%                               X.q = s
% cannot be solved as a linear system. 
% A singular value decomposition is used in the pinv algorithm to solve
% this type of problem. 
%
% Variables :
% The unknown vector q is composed of L elements.
% The data vector s is composed of M elements.
% The matrix X is composed of [LxM] elements. 
%
% If the problem is well-posed (M=L), the function is using a classical
% linear system solver.
% 
% Structure pp :
%
% Default :     pp = struct('regularization',0,...
%                           'doplot',1,...
%                           'compute_condition',1);
%
% pp.regularization : if 1, Tychonov regularization is made to improve the
%                     solution
% pp.doplot : if 1, plot the obtained data
% pp.compute_condition : if 1, the function returns the condition which is
% the ratio between the maximal sigular value and the minimal. This gives
% information on the quality of the pseudo inversion. 
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
    error('One matrix and one vector should be defined to solve the problem !')
end


nout = nargout - 1;
if nout > 0, pp.compute_condition = 1; end

varargout{1} = 0;

[M,tst] = size(s);

if tst > M
    s = s';
    [M,tst] = size(s);
end

if tst > 1,error('Data vector s must be 1-D array');end

[M_tst,L] = size(X);

if M_tst > L
    X = X';
    [M_tst,L] = size(X);
end

% check the size of the matrix compared to the data vector
if M_tst ~= M, error('Ymn : number of lines should be the same as p_vec'); end

% If M==L -> linear system
if M == L
    q = linsolve(X,s);

% Else -> ill posed problem    
else
    %%% Compute
    if pp.regularization 

        X_prod = X * (conj(X))';
        I = eye(size(X_prod));

        X = X_prod + alpha * I;

        if pp.compute_condition
            [~,Singular,~] = svd(X);
            max_S = max(max(Singular));
            min_S = min(min(Singular));
            condition = max_S/min_S;
            varargout{1} = condition;
        end

        Ymn_reg_inv = pinv(X);

        X_inv = (conj(X))' * Ymn_reg_inv ;

    else

        if pp.compute_condition
            [~,Singular,~] = svd(X);
            max_S = max(diag(Singular));
            min_S = min(diag(Singular));
            condition = max_S/min_S;
            varargout{1} = condition;
        end

        X_inv = pinv(X);    

    end

    q = X_inv * s;    

end

%%% plot
if pp.doplot
   h = figure('Name','Spherical Fourier Coefficients');
   plot(q,'k.','markersize',24)
   grid on
   title('Fourier coefficients $P_{mn}$','interpreter','latex')   
end

end