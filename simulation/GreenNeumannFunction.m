function Gn = GreenNeumannFunction(Rs,Rm,pp)
    
k = (2*pi*pp.freq) /pp.c;
SH_prod = 0;
Gn_tmp = 0;

r_m = Rm(1);
theta_m = Rm(2);
phi_m = Rm(3);

r_s = Rs(1);
theta_s = Rs(2);
phi_s = Rs(3);

for n = 0 : pp.maxOrder
    for m = -n : n
        mic_SH = getSphericalHarmonics(theta_m,phi_m,n,m);
        source_SH = getSphericalHarmonics(theta_s,phi_s,n,m);
        SH_prod = SH_prod + (mic_SH * conj(source_SH));
    end
    [hn_r0,~,~] = SphericalHankel1(n,k*r_s);
    [~,dhna,~] = SphericalHankel1(n,k*pp.a);
    
    Gn_tmp = Gn_tmp + ( (hn_r0/dhna) * SH_prod );
    SH_prod = 0;
end

Gn = -1/( k * (pp.a)^2 ) * Gn_tmp;

end