function plotit_XYZ(X, Y, Z, Xm, Ym, Zm, imin, Re, Rm)
%==========================================================================
% plotit_XYZ: Visualizes the 3D trajectories of the spacecraft and the moon 
%             in the Earth-Centered Inertial (ECI) frame, including an 
%             animation of their motion.
%
% Inputs:
%   X, Y, Z     - Spacecraft position vectors in ECI frame [km]
%   Xm, Ym, Zm  - Moon position vectors in ECI frame [km]
%   imin        - Index of the point of closest approach (perilune)
%   Re          - Radius of Earth [km]
%   Rm          - Radius of Moon [km]
%
% This function draws:
%   - Earth
%   - Moon at TLI, perilune, and end of simulation
%   - Spacecraft at TLI, perilune, and final position
%   - Animated trajectories of both spacecraft and Moon
%==========================================================================

    % Create figure window
    figure('Name','Trajectories of Spacecraft (red) and Moon (green)', ...
           'Color', 'w');
    
    % Generate a unit sphere mesh for plotting Earth and Moon
    [xx, yy, zz] = sphere(128);
    hold on
    
    % Plot coordinate axes in ECI frame (X, Y, Z)
    L = 20 * Re;  % Axis length scaling
    line([0 L], [0 0], [0 0], 'color','b');  % X-axis (blue)
    text(L, 0, 0, 'X', 'FontSize', 12, 'FontAngle', 'italic', 'FontName', 'Palatino');
    line([0 0], [0 L], [0 0], 'color','b');  % Y-axis (blue)
    text(0, L, 0, 'Y', 'FontSize', 12, 'FontAngle', 'italic', 'FontName', 'Palatino');
    line([0 0], [0 0], [0 L], 'color','b');  % Z-axis (blue)
    text(0, 0, L, 'Z', 'FontSize', 12, 'FontAngle', 'italic', 'FontName', 'Palatino');
    
    % Plot Earth as a semi-transparent sphere
    Earth = surfl(Re*xx, Re*yy, Re*zz);
    set(Earth, 'FaceAlpha', 0.5);
    shading interp
    
    % Plot spacecraft at TLI (departure point)
    plot3(X(1), Y(1), Z(1), 'o', 'MarkerEdgeColor', 'k', ...
        'MarkerFaceColor', 'k', 'MarkerSize', 5);
    
    % Plot spacecraft at perilune (closest approach to Moon)
    plot3(X(imin), Y(imin), Z(imin), 'o', 'MarkerEdgeColor', 'k', ...
        'MarkerFaceColor', 'k', 'MarkerSize', 4);
    
    % Plot spacecraft at end of simulation
    plot3(X(end), Y(end), Z(end), 'o', 'MarkerEdgeColor', 'r', ...
        'MarkerFaceColor', 'r', 'MarkerSize', 5);
    
    % Plot Moon at TLI position
    text(Xm(1), Ym(1), Zm(1), 'Moon at TLI');
    Moon1 = surfl(Rm*xx + Xm(1), Rm*yy + Ym(1), Rm*zz + Zm(1));
    set(Moon1, 'FaceAlpha', 0.99);
    shading interp
    
    % Plot Moon at perilune position
    Moon2 = surfl(Rm*xx + Xm(imin), Rm*yy + Ym(imin), Rm*zz + Zm(imin));
    set(Moon2, 'FaceAlpha', 0.99);
    shading interp
    
    % Plot Moon at end of simulation
    Moon3 = surfl(Rm*xx + Xm(end), Rm*yy + Ym(end), Rm*zz + Zm(end));
    set(Moon3, 'FaceAlpha', 0.99);
    shading interp
    
    % Initialize animated lines for spacecraft and Moon trajectories
    spacecraftLine = animatedline('Color','r','LineWidth',1.5);
    moonLine = animatedline('Color','g','LineWidth',1.2);
    
    % Main animation loop: plot the spacecraft and Moon moving over time
    for k = 1:length(X)
        % Update spacecraft trajectory
        addpoints(spacecraftLine, X(k), Y(k), Z(k));
        
        % Update Moon trajectory
        addpoints(moonLine, Xm(k), Ym(k), Zm(k));
        
        % Current position marker for spacecraft (triangle)
        h1 = plot3(X(k), Y(k), Z(k), '^', 'MarkerFaceColor', 'r', ...
            'MarkerEdgeColor', 'k', 'MarkerSize', 8);
        
        % Label for spacecraft
        t1 = text(X(k), Y(k), Z(k) + 0.1*Re, 'Spacecraft', ...
            'FontSize', 10, 'FontWeight','bold', ...
            'Color', 'r', 'HorizontalAlignment','center');
        
        % Current position marker for Moon (circle)
        h2 = plot3(Xm(k), Ym(k), Zm(k), 'o', 'MarkerFaceColor', 'g', ...
            'MarkerEdgeColor', 'k', 'MarkerSize', 8);
        
        % Label for Moon
        t2 = text(Xm(k), Ym(k), Zm(k) + 0.1*Re, 'Moon', ...
            'FontSize', 10, 'FontWeight','bold', ...
            'Color', 'g', 'HorizontalAlignment','center');
        
        % Update plot efficiently
        drawnow limitrate
        
        % Remove current markers and labels (to avoid clutter)
        if k < length(X)
            delete(h1); delete(h2);
            delete(t1); delete(t2);
        end
    end
end