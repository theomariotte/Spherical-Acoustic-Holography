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

% d�finition
global matlab_home;
matlab_home = 'H:\Mes documents\5A\Projet_5A\_GIT\Spherical-Acoustic-Holography\';

if ~exist(matlab_home,'dir')
   error('Initialization failed : MATLAB HOME does not exit'); 
end

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

% test existence des chemins
dd = true;
dd = dd && exist(data_dir,'dir') ...
        && exist(plot_dir,'dir')...
        && exist(simu_dir,'dir')...
        && exist(tool_dir,'dir')...
        && exist(test_dir,'dir');    
    
if dd == 0
   error('initialization failed : one path not found.'); 
end

% Ajout au PATH
addpath(data_dir);
addpath(plot_dir);
addpath(simu_dir);
addpath(tool_dir);
addpath(test_dir);


fprintf('Initialization done !\n')








