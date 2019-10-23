function [hn,diff_hn,r] = SphericalHankel1(n,r_in)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% [hn,diff_hn] = SphericalHankel1(n,r)
% 
% Function calculating the spherical Hankel of the first kind at the order
% N. It requires the calculation of the spherical Bessel functions of the
% first and second kinds. 
% The Hankel function is also differentiated using finite differences
% estimation of each Bessel function.
%
% see also SphericalBessel gradient
%
% Théo Mariotte - 16/10/2019 - S-NAH (ENSIM)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   
% spherical bessel function of the first and second kind at the order n
[jn,diff_jn,yn,diff_yn,r] = SphericalBessel(n,r_in);

% spherical Hankel function od the first kind defined as : hn(1) = jn + iyn
hn = jn + 1i*yn;

% derivative of the spherical hankel function 
diff_hn = diff_jn + 1j*diff_yn;


end