function Gn = getGreenNeumannFRF(ES_loc,mic_loc,Nmax,k,a)
% Compute the Green Neumann function for each microphone/source couple

% number of equivalent sources
nb_src = size(ES_loc,1);
% number of microphone
nb_mic = size(mic_loc,1);

% Build spherical harmonics and Hankel functions matrix for equivalent
% sources
idx = 1;
Y_mat_src = zeros(nb_src,(Nmax+1)^2);
hankel_src = zeros(nb_src,(Nmax+1)^2);

for n = 1 : Nmax
    [hn_src,~,~] = SphericalHankel2(n,k*ES_loc(:,1));
    [~,dhn_a,~] = SphericalHankel2(n,k*a);
   for m = -n:n
       Y_mat_src(:,idx) = getSphericalHarmonics(ES_loc(:,2),ES_loc(:,3),n,m);
       hankel_src(:,idx) = hn_src./dhn_a;
       idx = idx+1;    
   end
end

% Build spherical harmonics matrix for microphones
idx = 1;
Y_mat_mic = zeros(nb_mic,(Nmax+1)^2);
for n = 1 : Nmax
   for m = -n:n
       Y_mat_mic(:,idx) = getSphericalHarmonics(mic_loc(:,2),mic_loc(:,3),n,m);
       idx = idx+1;    
   end
end

% Compute the green neumann function for each microphone/source group
Gn = -1/(k*a^2) * (Y_mat_mic * (hankel_src .* Y_mat_src).');

end

