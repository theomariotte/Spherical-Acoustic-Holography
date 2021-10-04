# Spherical acoustic holography

Algorithmes permettant de reconstruire le champ acoustique à l'aide d'une antenne sphérique de microphones. 
Le champ de pression est mesuré sur une sphère rigide de rayon a. Les méthodes permettent de reconstruire la pression sur une 
sphère de rayon _r > a_. 

This repository proposes an implementation of spherical equivalent sources method [^1]. This method allows soundfield reconstruction outside of a rigid spherical microphone array. 

# Repository organization 

+ `init_spherical_NAH` : initializes Matlab pathes
+ `_tests` : some tests to validate implementations
+ `data` : data to be used for validation
+ `plot` : plotting functions
+ `simulation` : functions to generate simjulations of soundfield measurements
+ `tools` : tools for signal processing

[^1]: [1]E. Fernandez-Grande, « Sound field reconstruction using a spherical microphone array », The Journal of the Acoustical Society of America, vol. 139, nᵒ 3, p. 1168‑1178, mars 2016, doi: 10.1121/1.4943545.
