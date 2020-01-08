function handle = pressureMeasurementVisu(P_acou,Rm,Rs,pp_interp,pp_plot)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% handle = pressureMeasurementVisu(P_acou,Rm,Rs,pp_interp,pp_plot)
%
% Function allowing to visualize the sound field measured on a sphere by a
% set of microphones. The pressure array should be 1D with one pressure
% measurement per microphone. A field interpolation is performed before
% plotting the sound field in 3 different vusualizations :
%   - 'real' : the real part of the sound field is plotted
%   - 'abs'  : the norm of the sound field is plotted
%   - 'dB'   : the sound pressure level in dB SPL is computed
%
% Input :
%   - P_acou : array containing the sound pressure measurements
%   - Rm     : array containing microphones locations
%   - Rs     : array containing source(s) location(s)
%   - pp_interp : trcture containing interpolation parameters (radius of
%   the array, number of angles in one direction (along theta or phi),
%   interpolation method and plot or not the interpolation)
%   - pp_plot : plot parameters (see default values)
%
% Output :
%   - handle : handle of the figure containing the sound field
%   representation
%
% see also fieldInterpolation.m surf.m plot3.m
%
% Théo Mariotte - 11/2019 - ENSIM (SNAH project)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[P_acou,~] = shiftdim(P_acou);

if nargin < 5
    % default interpolation structure
    pp_plot = struct('showSource',1,...
                    'showMic',0,...
                    'plt_typ','dB',...
                    'clim',[40 100],...
                    'fontSize',12,...
                    'freq',[]);  
    warning('Default plot parameters used !');                
end

if nargin < 4
    % default interpolation structure
    pp_interp = struct('sphereRadius',15e-2,...
                'numAngle',30,...
                'interpType','natural',...
                'doplot',1);
    warning('Default interpolation parameters used !');
end

if nargin < 3, error('Not enough input argument !'); end

% field interpolation over a full grid 
[P_interp,x_grid,y_grid,z_grid] = fieldInterpolation(P_acou,Rm,pp_interp);

% compute the good sound pressure representation
switch pp_plot.plt_typ
    case 'real'
        P2plot = real(P_interp);
        if ~isempty(pp_plot.clim)
            clim = pp_plot.clim;
        else
            min_p = min(min(P_interp));
            max_p = max(max(P_interp));
            abs_lim = max(abs([min_p max_p]));
            clim = [-abs_lim abs_lim];
        end
        colorbarLabel = 'Sound Pressure (real) [Pa]';
        
    case 'dB'

        pref = 20e-6;
        P2plot = 20*log10(abs(P_interp)/pref);
        if ~isempty(pp_plot.clim)
            clim = pp_plot.clim;
        else
            clim = [min(min(P2plot)) max(max(P2plot))];
        end
        colorbarLabel = 'Sound Pressure Level [dB_{SPL}]';
    
    case 'abs'
        P2plot = abs(P_interp);       
        if ~isempty(pp_plot.clim)
            clim = pp_plot.clim;
        else
            clim = [min(min(P2plot)) max(max(P2plot))];
        end
        colorbarLabel = 'Sound Pressure (abs) [Pa]';
        
    otherwise 
        error('No plot type found - try (real) , (abs) or (dB)') 
end

% open the figure and plot
if nargout > 0
    handle = figure('Name','Sound pressure on the sphere');
else 
    handle = [];
end

hold on

% surface plot of the sound field
surf(x_grid,y_grid,z_grid,P2plot)      
colormap('jet')
shading('interp')
cc = colorbar;
cc.Label.String = colorbarLabel;
set(gca,'clim',clim);

% source plot
if pp_plot.showSource
    xs = Rs(:,1); ys = Rs(:,2); zs = Rs(:,3);
    plot3(xs,ys,zs,'ko','linewidth',3)
end

% transducers plot
if pp_plot.showMic
    x = Rm(:,1); y = Rm(:,2) ; z = Rm(:,3);
    plot3(x,y,z,...
        'linestyle','none',...
        'Marker','.',...
        'linewidth',5,...
        'color',[.5 .5 .5],...
        'markersize',12);
end

hold off

% label titles
xlabel('X','interpreter','latex')
ylabel('Y','interpreter','latex')
zlabel('Z','interpreter','latex')
if isempty(pp_plot.freq)
    title(sprintf('Measured sound pressure (%s)',pp_plot.plt_typ),...
        'interpreter','latex')
else
    title(sprintf('Measured sound pressure (%s) : f = %d Hz',pp_plot.plt_typ,pp_plot.freq),...
        'interpreter','latex')    
end

axis('equal')
grid on


