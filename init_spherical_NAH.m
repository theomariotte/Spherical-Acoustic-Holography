%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Initialisation de la méthode d'holographie acoustique sphérique.
% Dossiers initialisés :
%   - spherical functions 
%   - plot
%   - tools (résolution problème mal posé)
%   - data (positions micros, mesures)
%   - tests (tous les mains pour tester les différentes parties)
%   - simulation (simulation pour tester la méthode)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear; clc; 
close all;

% définition
global matlab_home;
matlab_home = 'C:\Users\Théo\Documents\1_WORK\01_ENSIM\5A\Projet 5A\Spherical-Aoustic-Holography\';

% Chemins en temps que variables globales
global data_dir;
global plot_dir;
global simu_dir;
global tool_dir;
global test_dir;

data_dir = [matlab_home 'data\'];
plot_dir = [matlab_home 'plot\'];
simu_dir = [matlab_home 'simulation\'];
tool_dir = [matlab_home 'tools\'];
test_dir = [matlab_home '_tests\'];

% Ajout au PATH
dd = true;

dd = dd && exist(data_dir,'dir') ...
        && exist(plot_dir,'dir')...
        && exist(simu_dir,'dir')...
        && exist(tool_dir,'dir')...
        && exist(test_dir,'dir');    
    
if dd == 0
   error('initialization failed : one path not found.'); 
end

addpath(data_dir);    
addpath(plot_dir);
addpath(simu_dir);
addpath(tool_dir);
addpath(test_dir);


fprintf('Initialization done !\n')








