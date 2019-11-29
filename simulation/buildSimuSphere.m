function Rm = buildSimuSphere(ntheta,nphi,a,doplot)

theta_grid = linspace(-pi/2,pi/2,ntheta)';
phi_grid = linspace(-pi,pi,nphi)';
k = 1;
for ii = 1 : nphi
    for jj = 1 : ntheta
        x(k) = a * cos(theta_grid(jj)) .* cos(phi_grid(ii));
        y(k) = a * cos(theta_grid(jj)) .* sin(phi_grid(ii));
        z(k) = a * sin(theta_grid(jj));
        k = k + 1;
    end
end

Rm = [x' y' z'];

if doplot
figure
plot3(x,y,z,'k.','linewidth',3,'markersize',16)
axis equal
title(sprintf('Antenne avec M = %d microphones',length(x)))
grid on
xlabel('X')
xlabel('Y')
xlabel('Z')
end
end