function G = GreenFunction(r,r0,freq,pp_source)

if nargin < 4
    warning('No source parameters defined. Unitary solution only.')
    pp_source = struct('Q',1e-3,...
                       'c',340,...
                       'rho',1,...
                       'unitary',1);                     
end

if nargin < 3, error('Not enough parameters.'); end
% Wave number at the current frequency
w = 2*pi*freq;
k = w / pp_source.c;

% source-receiver distance
R = norm(r - r0);
G = exp(1i*k*R)/(4*pi*R);

% Normalisation of the Green's function in case of a monopole source
if ~pp_source.unitary
    G = 1j * k * pp_source.rho * pp_source.c * pp_source.Q * G;
end

end