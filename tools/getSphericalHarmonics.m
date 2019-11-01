function [Y] = getSphericalHarmonics(theta,phi,n,m,pp)  

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% [Y] = getSphericalHarmonics(theta,phi,n,m,pp)  
% 
% Function calculating the spherical harmonic for given degree N and order
% M.
% Theta is the azimuth angle. Phi is the elevation angle. 
% The structure pp contains different parameters which are defined below :
%   pp = struct('norm',1,... 
%               'doplot',1);
% norm : normalisation de l'harmonique sphérique
% doplot : plot the associated legendre function (elevation function)
%
% see also legendre 
%
% Théo Mariotte - 15/10/2019 - S-NAH (ENSIM)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Function initialization
if nargin < 5
    pp = struct('norm',1,...
                'doplot',0);
end

if nargin < 4
    warning('Function order not defined. Default value m = 0 used')
    m = 0;
end

if nargin < 3
    warning('Function degree not defined. Default value n = 0 used')
    n = 0; 
end

if nargin < 2, error('Both elevation - azimuthal angles grid should be defined.'); end
if size(theta) ~= size(phi), error('Theta and Phi grids should be the same size'); end
if abs(m) > n, error('the order condition |m| <= n should be verified.');end
   
% check if m is negative and odd
isneg = m < 0;
isodd = mod(m,2) == 1;

% resize each grid as a column array
sz = size(theta);
theta = theta(:);
phi = phi(:);

% set m as a positive value
if isneg, m = abs(m); end

% calculation of the legendre associated functions for the degree n. All
% the orders m = 0,...,n are calculated. 
Pmn_whole = legendre(n,cos(theta));

% extract only the useful part (m order) for spherical harmonic calculation
Pmn = Pmn_whole(m+1,:)';

% Spherical harmonic normalization to respect the orthogonality principle
if pp.norm
    num = (2*n+1) * factorial(n-m);
    den = 4*pi * factorial(n+m);
    C = sqrt(num/den);
else
    C = 1; 
end

% current spherical harmonic (complex)
Y = C * Pmn .* exp(-1j*m*phi);

% in case of negative m, the spherical harmonic can be computed from the
% one at the same m but positive (see Williams Eq (6.44) pp 191)
if isneg
    Y = conj(Y); 
    if isodd
        Y = -Y;
    end
end

% reshape at the same size as the original grid
Y = reshape(Y,sz);

% plot the current associated legendre function
if pp.doplot
    hh = figure('Name',sprintf('Associated Legendre function : degree %d ; order %d',n,m));
    plot(theta(1:sz(1)),Pmn(1:sz(1)),'k','linewidth',2)
    xlabel('Elevation angle $\theta$','interpreter','latex')
    title(sprintf('Associated Legendre function $P_{%d}^{%d}$',n,m),'interpreter','latex')
    grid on
end


end