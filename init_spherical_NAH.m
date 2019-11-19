%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Initialisation de la m�thode d'holographie acoustique sph�rique.
% Dossiers initialis�s :
%   - spherical functions 
%   - plot
%   - tools (r�solution probl�me mal pos�)
%   - data (positions micros, mesures)
%   - tests (tous les mains pour tester les diff�rentes parties)
%   - simulation (simulation pour tester la m�thode)
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

% d�finition
matlab_home = 'H:\Mes documents\5A\Projet_5A\_GIT\Spherical-Acoustic-Holography\';

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








