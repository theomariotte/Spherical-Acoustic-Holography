function [Xeq,Yeq,Zeq] = equivalentSourcesGrid(R,Ns,xlim,ylim,zlim)

if nargin < 4    
    rand_grid  = rand(Ns);
    phi_grid = rand_grid * 2*pi;
    theta_grid = rand_grid * pi;
    radius_grid = R * ones(size(rand_grid));        
    
    [Xeq,Yeq,Zeq] = cartesianCoordinates(radius_grid,phi_grid-pi,theta_grid);
    
elseif nargin < 2
    
else
    
end



end