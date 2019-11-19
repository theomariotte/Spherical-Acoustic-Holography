9function RadialFuncVisu(hn,r,pp)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% RadialFuncVisu(hn,diff_hn,r,pp)
%
% Function to visualize spherical radial functions (Hankel of the first
% kind, Bessel of first and second kind). the derivative can aslo be
% visualize if specified in pp. 
%
% Théo Mariotte - 17/10/2019
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


if nargin < 3
   pp = struct('diff_hn',[],'order',NaN,'fontsize',12,...
               'x_label','Radius r'); 
end

doplotdiff = ~isempty(pp.diff_hn);

%%%%%%%%%%%%%%%%%%%%%%%%%% Spherical hankel %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

h1 = figure('Name','First kind spherical Hankel');
% create axes
ax1 = subplot(2,1,1);
ax2 = subplot(2,1,2);
% holding on
hold(ax1,'on');
hold(ax2,'on');

% plot on each axis
plot(ax1,r,real(hn),'k','linewidth',2)
plot(ax2,r,imag(hn),'k','linewidth',2)

% if the derivative has to be plotted
if doplotdiff
    dhn = pp.diff_hn;
    plot(ax1,r,real(dhn),'k','linewidth',2,'linestyle','--')
    plot(ax2,r,imag(dhn),'k','linewidth',2,'linestyle','--')
      
    % legend (only in this case because more than one curve per axe)
    legend(ax1,{'$h_n$','$d h_n$'},'interpreter','latex')
    legend(ax2,{'$h_n$','$d h_n$'},'interpreter','latex')
end

% set some parameters
set(ax1,'fontsize',pp.fontsize)  
set(ax2,'ylim',[-50 50],'fontsize',pp.fontsize)  
% set grids on
grid(ax1,'on'); 
grid(ax2,'on');
% X labels
xlabel(ax1,pp.x_label,'interpreter','latex');
xlabel(ax2,pp.x_label,'interpreter','latex');
% Y labels
ylabel(ax1,'$Re(h_n) = j_n$','interpreter','latex'); 
ylabel(ax2,'$Im(h_n) = y_n$','interpreter','latex');
% Titles
title(ax1,sprintf('Spherical Hankel Function (real) : $n=%d$',pp.order),'interpreter','latex')
title(ax2,sprintf('Spherical Hankel Function (imag) : $n=%d$ : ',pp.order),'interpreter','latex')
% hold off
hold(ax1,'off');
hold(ax2,'off');



end