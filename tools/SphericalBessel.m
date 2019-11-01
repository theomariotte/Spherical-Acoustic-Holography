function [jn,diff_jn,yn,diff_yn,r] = SphericalBessel(n,r_in)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% [jn,diff_jn] = SphericalBesselj(n,r)
% 
% Function calculating the spherical Bessel functions from order 0 to Nmax. 
% The differential along the r axis is also estimated using finite 
% differences.
%
% see also besselj bessely gradient
%
% Théo Mariotte - 16/10/2019 - S-NAH (ENSIM)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if nargin < 2, error('Rdius array should be defined');end

% keep values of r strictly higher than 0
r = r_in(r_in > 0);
if isempty(r), error('r should be an array of positive values');end

% distance between points along r axis
dr = mean(diff(r));
    
%%% current first kind Bessel function
J = besselj(n+1/2,r);    

% Spherical Bessel function
jn = sqrt(pi./(2*r)) .* J;

%%% current second kind Bessel function
Y = bessely(n+1/2,r);    

% Spherical Bessel function
yn = sqrt(pi./(2*r)) .* Y;

% derivatives
% Derivatives of bessel functions are computed using recursive relations
% between functions of different orders. Only the order n=0 is computed
% using the analytical expression. 
if n == 0
    % analytical derivatives for the 0 order
    diff_jn = (r .* cos(r) - sin(r))./(r.^2);
    diff_yn = (r .* sin(r) + cos(r))./(r.^2);
else
    % recursive relations
    [jn_prev,~,~,~,~] = SphericalBessel(n-1,r);
    [~,~,yn_prev,~,~] = SphericalBessel(n-1,r);
    
    diff_jn = jn_prev - (n+1)./r .* jn;
    diff_yn = yn_prev - (n+1)./r .* yn;
    
end


end