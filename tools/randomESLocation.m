function [Xs,Ys,Zs] = randomESLocation(R,Ns,plotfl)
% [Xs,Ys,Zs] = randomESLocation(R,Ns)
% Function that computes the coordinates of equivalent sources over a
% sphere. This is used when actual source location is unknown. 
% Coordinates are computed using a uniform law.
%
% Input :
%   - R : radius of the sphere where equivalent sources are located
%   - Ns : number of sources
%   - plotfl : plot flag (>1 if you want to plot)
%
% Output :
%   - [Xs,Ys,Zs] : sources coordinates in cartesian coordinates


if nargin < 3, plotfl=0;end

theta = randn(Ns,1)*pi;
phi = randn(Ns,1)*2*pi;
rho = ones(Ns,1) * R;

[Xs,Ys,Zs] = cartesianCoordinates(rho,theta,phi);

if plotfl
   h = figure('Name','Random equivalent sources locations');
   plot3(Xs,Ys,Zs,...
       'linestyle','none',...
       'marker','o',...
       'linewidth',2,...
       'color','k');
   xlabel('X')
   ylabel('Y')
   zlabel('Z')
   axis equal
end

end