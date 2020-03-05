function Pmn = sphericalTFSolver(Y,pt)

nmic = size(Y,1);
nb_SH = size(Y,2);

if nmic < nb_SH
    error('the problem will not be solved correctly if sz(Y,1) < sz(Y,2). order is too big !')
end

Y_inv = pinv(Y);

Pmn = Y_inv * pt;

end