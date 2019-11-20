function flag = printFigFmt(fig_handle,path,fname,ext)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% flag = printFigFmt(fig_handle,path,fname,ext)
%
% Foction permettant d'enregistrer une image selon diff�rents formats
% ('png','eps','emf').
% 
% Entr�es : 
%   - fig_handle : handle de la figure � sauvegarder
%   - path       : chemin du dossier o� sauvegarder l'image
%   - fname      : nom du fichier � enregistrer. S'il existe d�j�, une
%   nouvelle version avec '_0x' au bout sera ajout�e.
%   - ext        : extension du fichier sous lequel il doit �tre enregistr�
%
% Sorties :
%   - flag       : bool�en informant sur le succ�s ou non de
%   l'enregistrement
%
% T. MARIOTTE - 20/11/2019 - ENSIM
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
try 
    flag = 1;
    outImName = getOutFileName(path,fname,ext);
    
    if strcmp(ext,'eps') == 1 || strcmp(ext,'.eps') == 1
        fmt = '-depsc';
    elseif strcmp(ext,'png') == 1 || strcmp(ext,'.png') == 1
        fmt = '-dpng';
    elseif strcmp(ext,'emf') == 1 || strcmp(ext,'.emf') == 1
        fmt = '-dmeta';
    else
        fprintf('No image format found !\n');
        flag = 0;
    end
    
    print(fig_handle,outImName,fmt);
catch
   fprintf('Image saving failed !\n')
   flag = 0;
end

end