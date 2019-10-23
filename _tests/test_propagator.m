% test propagator

clear; clc; 
close all


a = 15e-2;
r = 2*a;
k = linspace(0.001,100,2048);

n_vec = [4 7 9];

h = figure('Name','Propagator');
hold on
for i = 1 : length(n_vec)
    
    n = n_vec(i);
    
    % derivative of the spherical hankel function on the surface of the
    % sphere
    [~,dhn_a,~] = SphericalHankel1(n,k*a);

    % first and second kind bessel function at the microphone location
    [jn_r,~,yn_r,~,~] = SphericalBessel(n,k*r);

    % derivatives of Bessel functions on the surface of the sphere
    [~,djn_a,~,dyn_a,~] = SphericalBessel(n,k*a);

    % propagator
    Gn = (k*a).^2 .* (jn_r .* dyn_a - djn_a .* yn_r);
    
    logGn = 20*log10(abs(Gn));
    plot(k*a,logGn,'linewidth',2);    
    L{i} = sprintf('n = %d',n);
end
set(gca,'ylim',[-20 80])
grid on
xlabel('ka')
legend(L);