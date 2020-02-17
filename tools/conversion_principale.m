type_data=0; % Pressions acoustiques

no_devices=1:3; % n°carte
no_voies{1}=1:16; % n°voie
no_voies{2}=1:16;
no_voies{3}=1:4;
nbvoies=36;


chem='D:\Users\pxi\Bureau\Theo';
addpath(chem)

% Le fichier .mat existe
% Attention dans ce cas, se mettre dans le répertoire du fichier
% [DATA,t,fe]=convert_data_TDMS(nbvoies,type_data,no_devices,no_voies,0);

% Le fichier .mat n'existe pas
[DATA,t,fe]=convert_data_TDMS(nbvoies,type_data,no_devices,no_voies,1);

visu_data_TDMS(t,DATA,nbvoies,1)
