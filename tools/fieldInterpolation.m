function [P_interp,x_grid,y_grid,z_grid] = fieldInterpolation(P,mic_pos,param)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
% P_interp = fieldInterpolation(P,mic_pos,numAngle)
%
% Function to interpolated a given field on a sphere to get a better
% resolution for plotting. 
%
% Inputs : 
%   * P = array containing the field to be interpolated. 
%   * mic_pos = tranducers locations arranged in a [Mx3] array
%   * param = structure containing interpolation parameters
%
% Ouputs :
%   * P_interp = interpolated field on the new grid
%   * x_grid = new grid on the X axis
%   * y_grid = new grid on the Y axis
%   * z_grid = new grid on the Z axis
%
% Théo Mariotte - 11/2019 - ENSIM (SNAH project)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % microphone's locations in cartesian coordinates        
    x = mic_pos(:,1);
    y = mic_pos(:,2);
    z = mic_pos(:,3);
    
    % Interpolate the sound field
    F = scatteredInterpolant(x,y,z,P,param.interpType);
    
    % Build the interpolation grid
    nTheta0 = param.numAngle; 
    nPhi0 = nTheta0;
    az_ = linspace(-pi, pi, nPhi0);
    elev_   = linspace(-pi/2, pi/2, nTheta0);
    [az_q, elev_q] = meshgrid(az_, elev_);
    
    % Constant radius sphere
    radius = param.sphereRadius*ones(nTheta0, nPhi0);
    
    % Interpolation grid converted in cartesian coordinates
    [x_grid,y_grid,z_grid] = sph2cart(az_q, elev_q, radius);
    
    % Interpolated pressure
    P_interp = F(x_grid,y_grid,z_grid);
    
    
    if param.doplot
        
        % plot to evaluate the grid
%         figure('Name','Interpolation grid')
%         hold on;
%         surf(x_grid,y_grid,z_grid,ones(size(x_grid)));
%         colormap('gray');
%         set(gca,'clim',[0 1]);
%         plot3(x,y,z,'rx','linewidth',3);
%         hold off;
%         axis equal;
%         xlabel('X');
%         ylabel('Y');
%         zlabel('Z');
%         title('Field interpolation : grid evaluation');
%         legend('Interpolation grid','Transducers');
        
        % Compare initial sound pressuer with the interpolated one
        [~,phi_init,theta_init] = sphericalCoordinates(x,y,z);        
        [~,az_q,elev_q] = sphericalCoordinates(x_grid,y_grid,z_grid);        
        figure('Name','Sound FieldInterpolation');
        hold on
        surf(az_q,fliplr(elev_q),abs(P_interp))
        c = colorbar;
        c.Label.String = 'Interp. sound pressure (abs) [Pa]';
        shading('Interp')
        plot3(phi_init,...
              theta_init,...
              abs(P),...
              'k.',...
              'linewidth',4);
          
          % axis labels
        xlabel('Azimuth $\phi$','interpreter','latex')
        ylabel('Elevation $\theta$','interpreter','latex')
        zlabel('Sound pressure $|p(a,\theta,\phi)|$','interpreter','latex')
                       
        grid on
        legend('Interpolated','Measurement');
        set(gca,'fontsize',12);
        view(80,50);
    end
  
end

