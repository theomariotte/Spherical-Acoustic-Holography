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


% Chemins en temps que variables globales
global matlab_home;
global data_dir;
global plot_dir;
global simu_dir;
global tool_dir;
global test_dir;

% définition
matlab_home = 'C:\Users\Théo\Documents\1_WORK\01_ENSIM\5A\Projet 5A\Spherical-Aoustic-Holography\';

data_dir = [matlab_home 'data\'];
plot_dir = [matlab_home 'plot\'];
simu_dir = [matlab_home 'simulation\'];
tool_dir = [matlab_home 'tools\'];
test_dir = [matlab_home '_tests\'];

% Ajout au PATH
try
    addpath(data_dir);
    addpath(plot_dir);
    addpath(simu_dir);
    addpath(tool_dir);
    addpath(test_dir);
catch 
    fprintf('Path initialization failed. Check matlab_home name !\n')
end

fprintf('Initialization done !\n')








