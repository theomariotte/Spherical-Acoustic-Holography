
Nmax = 10;
ka = 2;

f = 500;
k = 2*pi*f /343;

a = ka/k;
r0 = 5*a;
phi0 = 0;
theta0 = 0;

idx = 1;
for n = 0 : Nmax
    [hn_r0,~,~] = SphericalHankel1(n,k*r0);
    [~,dhn_a,~] = SphericalHankel1(n,ka);
    P0n(n+1) = -sqrt(4*pi*(2*n + 1)) * hn_r0 / (((ka)^2 * dhn_a));
    for m = -n : n
        Ymn = getSphericalHarmonics(theta0,phi0,n,m);
        Pmn(idx) = -4*pi*( hn_r0 / ((ka)^2 * dhn_a) ) * conj(Ymn);
        idx = idx + 1;
    end
end

nn= 0 : length(P0n)-1;
figure
hold on
plot(nn,20*log10(abs(P0n)),'-r','linewidth',3)
plot(nn,20*log10(abs(Pmn(Pmn ~=0))),'xb','linewidth',3)
grid on
set(gca,'xlim',[0 30],'ylim',[-80 0])
xlabel('n')
ylabel('20log_{10}(Pmn)')
title('Comparison of Fourier coefficients')

