function P_i = incidentPressure(r,r0,pp_simu)

c = pp_simu.freq;
f = pp_simu.freq;
w = 2*pi*f;
k = w/c;

rho = pp_simu.rho;
Q = pp_simu.Q;

R = norm(r - r0);

G = exp(1i*k*R)/(4*pi*R);

P_i = (-1i*rho*c*Q) * G;


end