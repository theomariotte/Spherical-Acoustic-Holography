function visu_data_TDMS(t,DATA,nbvoies,fig)
%function visu_data_TDMS(t,DATA,nbvoies,fig)
% Visualisation des données
% t : vecteur temps
% DATA : données (attention le nom DATA doit être employé)
% nbvoies : nbvoies
% fig : numéro de figure pour l'affichage
nbfig=floor(nbvoies/4);
nbfigbis=rem(nbvoies,4);
nomdata='DATA.voie';

for i=1:nbfig
    figure(fig+i-1)
    subplot(411)
    plot(t,eval([nomdata,num2str((i-1)*4+1)]))
    subplot(412)
    plot(t,eval([nomdata,num2str((i-1)*4+2)]))
    subplot(413)
    plot(t,eval([nomdata,num2str((i-1)*4+3)]))
    subplot(414)
    plot(t,eval([nomdata,num2str((i-1)*4+4)]))
    xlabel('Temps [s]')
    set(gcf,'name',['Voies ',num2str((i-1)*4+1),'-',num2str(i*4)])
end

if (nbfigbis~=0)
    figure(fig+nbfig)
    for i=1:nbfigbis
        feval('subplot',['41',num2str(i)])
        plot(t,eval([nomdata,num2str(nbfig*4+i)]))
    end
    xlabel('Temps [s]')
    if (nbfigbis==1)
        set(gcf,'name',['Voies ',num2str(nbfig*4+1)])
    else
        set(gcf,'name',['Voies ',num2str(nbfig*4+1),'-',num2str(nbfig*4+nbfigbis)])
    end
end