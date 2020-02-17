function [DATA,t,fe]=convert_data_TDMS(nbvoies,type_data,no_devices,no_voies,choix_conv,nomfich)
%function [DATA,t,fe]=convert_data_TDMS(nbvoies,type_data,no_devices,no_voies,choix_conv,nomfich)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Permet de r�cup�rer les donn�es dans la matrice DATA
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% nbvoies : nombre de voies (du m�me type)
% type_data : type de donn�es (pressions, acc�l�rations, tensions)
% type_data=0; % Pressions acoustiques
% type_data=1; % Acc�l�rations
% type_data=2; % Tensions
% no_devices : num�ros des cartes
% -> exemple : no_devices=[1:4];
% no_voies : num�ros des voies pour chaque carte
% -> exemple :
% no_voies{1}=1:16;
% no_voies{2}=1:16;
% no_voies{3}=1:16;
% no_voies{4}=1;
% choix_conv : choix de la conversion -> conversion totale : 1
%              -> conversion partielle : 0 (le fichier .mat existe)
% si choixconv=0 pas d'appel � simpleConvertTDMS ->
% Note : la fonction doit �tre lanc�e dans le r�pertoire des donn�es
% nomfich : nom d'un fichier pour la conversion
% DATA : cellule contenant les signaux (DATA.voie1, DATA.voie2...)
% t : vecteur temps
% fe : fr�quence d'�chantillonnage

nb_devices=length(no_devices);
device='PXI_4496_0';

switch type_data
    case 0
        nomcapteur='Pression_acoustique';
        nomchaine='Pressionacoustique';
    case 1
        nomcapteur='Accel';
        nomchaine='Accel';
    case 2
        nomcapteur='Tension';
        nomchaine='Tension';
end


% Chemin de programmes de conversion
chem='D:\Users\All Users\Documents\pb_conversion\humphreysb-ConvertTDMS-bf58f6b';

addpath(chem)

% Conversion d'un fichier TDMS
% par cr�ation d'un fichier .mat
% Ne fonctionne qu'avec des versions de Matlab > 7.04
% Tr�s bien avec Matlab 7.7 (R2008b)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if choix_conv
    if (nargin<6)
        matFileName=simpleConvertTDMS
    else
        matFileName=simpleConvertTDMS(nomfich)
    end
    % permet de travailler dans le r�pertoire du fichier acquis
    nomfich=eval('matFileName{1}');
else
    nomfich=[pwd,'\',nomcapteur,'.mat']
    end
%matFileName={'D:\Users\pxi\Bureau\Justine\data\16032017_155712\Pression_acoustique.mat'};

%nomfich2=nomfich(1:end-(length(nomcapteur)+4));
position=strfind(nomfich,'/');
nomfich2=nomfich(1:position(end));
feval('cd',nomfich2)

feval('load',nomcapteur)
% load Pression_acoustique
% load Accel

% R�cup�ration de la date sous la forme
%   date_mes='12012016150839' pour 12-Jan-2016 15:08:39
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Strat�gie 2 : aller chercher la date dans ci.Objectn.long_name
% de la forme 12/01/2016 15:08:39
it=0;
longueur=20;
date_mes=ci.Object1.long_name;
while ((length(date_mes)<longueur) & (it<36))
    it=it+1;
    eval(['date_mes=ci.Object',num2str(it+1),'.long_name;'])
end
date_mes=[date_mes(1:2) date_mes(4:5) date_mes(7:10) date_mes(12:13),...
    date_mes(15:16) date_mes(18:19)];
% Fin strat�gie 2

% R�cup�ration de la fr�quence d'�chantillonnage
fe=1/Root.Property.logdt;
disp(['Fr�quence d''�chantillonnage : ',num2str(fe),' Hz'])

% R�cup�ration des noms des variables (structures) contenant les donn�es
% dans la variable nom_struct
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

n=0;

for i=1:nb_devices
    nbvoies_dev=length(no_voies{i}); % nbvoies_dev : nombre de voies par device
    for j=1:nbvoies_dev
        n=n+1;
        voie{n}=[device,num2str(no_devices(i)),'_ENS_ai',num2str(no_voies{i}(j)-1)];
    end
end


disp('Noms des structures')
disp('*******************')
for i=1:nbvoies 
    nom_voie=voie{i};
    s=['d',date_mes,nomchaine,'AllData',nom_voie];
    ls=length(s);
    it=0;
    while (~exist(eval('s'))&(it<3))
        it=it+1;
        s(ls+1)=num2str(it);
    end
    nom_struct{i}=s;
    disp(nom_struct{i})    
end

% R�cup�ration des donn�es
don=eval(eval('s'));
nbech=don.Total_Samples;

for i=1:nbvoies
    don=eval(nom_struct{i});
    eval(['DATA.voie',num2str(i),'=don.Data;'])
    %DATA(i,:)=don.Data.';
    %feval('clear',nom_struct{i})
end

% Elaboration du vecteur temps
t=(0:nbech-1)/fe;

