function Gn = GreenNeumannFunction(Rs,Rm,freq,Nmax,pp_source)
    

k = (2*pi*freq) / pp_source.c;
SH_prod = 0;
Gn_tmp = 0;

r_m = Rm(1);
theta_m = Rm(2);
phi_m = Rm(3);

r_s = Rs(1);
theta_s = Rs(2);
phi_s = Rs(3);

for n = 0 : Nmax
    for m = -n : n
        mic_SH = getSphericalHarmonics(theta_m,phi_m,n,m);
        source_SH = getSphericalHarmonics(theta_s,phi_s,n,m);
        SH_prod = SH_prod + (mic_SH * conj(source_SH));
    end
    [mic_jn,~,~,~,~] = SphericalBessel(n,k*r_m);
    [~,source_diff_jn,~,~,~] = SphericalBessel(n,k*r_s);
    
    Gn_tmp = Gn_tmp + ( (mic_jn/source_diff_jn) * SH_prod );
end

Gn = 1/( k * r_s^2 ) * Gn_tmp;

end