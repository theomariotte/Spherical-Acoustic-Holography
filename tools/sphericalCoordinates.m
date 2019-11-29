function [rho,phi,theta] = sphericalCoordinates(x,y,z)
    % conversion sphérique
    rho = sqrt(x.^2 + y.^2 + z.^2);
    phi = atan2(y,x);
    theta = acos(z./rho);
end