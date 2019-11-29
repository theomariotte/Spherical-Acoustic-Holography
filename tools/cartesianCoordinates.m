function [x,y,z] = cartesianCoordinates(rho,phi,theta)

x = rho .* cos(theta) .* cos(phi);
y = rho .* cos(theta) .* sin(phi);
z = rho .* sin(theta);

end