% Reconstruction test from experimental data

clear; clc; 
close all;

%% load data measured by the array

d = 41.5e-2;
dataPath = ['/Users/theomariotte/Documents/01_work/ENSIM/5A/Projet/'...
    'Spherical-Acoustic-Holography/data/12022020_110356/'];
fname = 'data/12022020_110356/Pression_acoustique.mat';

type_data=0; % Pressions acoustiques

no_devices=1:3; % n°carte
no_voies{1}=1:16; % n°voie
no_voies{2}=1:16;
no_voies{3}=1:4;
nbvoies=36;

addpath(dataPath);

[DATA,t,fs]=convert_data_TDMS(nbvoies,type_data,no_devices,no_voies,0);

%% compute spectrum

for k = 1 : nbvoies
   
    crnt_x = eval(sprintf('DATA.voie%d',k));
    Gxx = pwelch(crnt_x,hann(M),M/2,M); 
    % comment obtenir un spectre complexe ?? 
    
end

%% compute equivalent sources locations



%% compute sources strenght

%% 
