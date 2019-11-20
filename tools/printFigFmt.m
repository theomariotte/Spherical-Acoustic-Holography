function flag = printFigFmt(fig_handle,path,fname,ext)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% flag = printFigFmt(fig_handle,path,fname,ext)
%
% Foction permettant d'enregistrer une image selon différents formats
% ('png','eps','emf').
% 
% Entrées : 
%   - fig_handle : handle de la figure à sauvegarder
%   - path       : chemin du dossier où sauvegarder l'image
%   - fname      : nom du fichier à enregistrer. S'il existe déjà, une
%   nouvelle version avec '_0x' au bout sera ajoutée.
%   - ext        : extension du fichier sous lequel il doit être enregistré
%
% Sorties :
%   - flag       : booléen informant sur le succès ou non de
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