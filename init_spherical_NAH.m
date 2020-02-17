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
% matlab_home = 'H:\Mes documents\5A\Projet_5A\_GIT\Spherical-Acoustic-Holography\';
% matlab_home = 'C:\Users\Th�o\Documents\1_WORK\01_ENSIM\5A\Projet 5A\Spherical-Aoustic-Holography\';
matlab_home = '/Users/theomariotte/Documents/01_work/ENSIM/5A/Projet/Spherical-Acoustic-Holography/';
if ~exist(matlab_home,'dir')
   error('Initialization failed : MATLAB HOME does not exit'); 
end

% Chemins en temps que variables globales
global data_dir;
global plot_dir;
global simu_dir;
global tool_dir;
global test_dir;

data_dir = [matlab_home 'data'];
plot_dir = [matlab_home 'plot'];
simu_dir = [matlab_home 'simulation'];
tool_dir = [matlab_home 'tools'];
test_dir = [matlab_home '_tests'];

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
fs=44100;
duration = 0.25;
N = ceil(fs*duration);
tt = (0:N-1)/fs;
x = chirp(tt,400,duration,1000);
x2 = chirp(tt,1000,duration,400);
soundsc([x x2],fs);






