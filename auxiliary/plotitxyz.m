function plotitxyz(x, y, z, xm, ym, zm, imin, Re, Rm)
%==========================================================================
% plotitxyz: Plots the trajectory of the spacecraft in a Moon-fixed 
%            rotating reference frame, along with the Moon and Earth.
%
% Inputs:
%   x, y, z      - Spacecraft position vectors in Moon-fixed frame [km]
%   xm, ym, zm   - Moon position vectors in Moon-fixed frame [km]
%   imin         - Index of closest approach (perilune)
%   Re           - Radius of Earth [km]
%   Rm           - Radius of Moon [km]
%
% The figure shows:
%   - Spacecraft trajectory (red line)
%   - Moon trajectory (green line)
%   - Earth (semi-transparent sphere)
%   - Moon at TLI, perilune, and final position
%   - Coordinate axes of the Moon-fixed frame
%==========================================================================

    % Create figure window
    figure('Name', 'Spacecraft trajectory in Moon-fixed rotating frame', ...
        'Color', [1 1 1]);

    % Generate unit sphere mesh for Earth and Moon rendering
    [xx, yy, zz] = sphere(128);
    hold on

    % Plot spacecraft trajectory (red curve)
    plot3(x, y, z, 'r', 'LineWidth', 2.0);

    % Plot Moon's trajectory (green curve)
    plot3(xm, ym, zm, 'g', 'LineWidth', 0.5);

    % Plot Earth as a semi-transparent sphere at the origin
    Earth = surfl(Re * xx, Re * yy, Re * zz);
    set(Earth, 'FaceAlpha', 0.5);
    shading interp

    % Draw coordinate axes of Moon-fixed frame
    L1 = 63 * Re;  % x-axis length scaling
    L2 = 20 * Re;  % y-axis length scaling
    L3 = 29 * Re;  % z-axis length scaling
    line([0 L1], [0 0], [0 0], 'color', 'k');  % x-axis (black)
        text(L1, 0, 0, 'x', 'FontSize', 12, 'FontAngle', 'italic', 'FontName', 'Palatino');
    line([0 0], [0 L2], [0 0], 'color', 'k');  % y-axis (black)
        text(0, L2, 0, 'y', 'FontSize', 12, 'FontAngle', 'italic', 'FontName', 'Palatino');
    line([0 0], [0 0], [0 L3], 'color', 'k');  % z-axis (black)
        text(0, 0, L3, 'z', 'FontSize', 12, 'FontAngle', 'italic', 'FontName', 'Palatino');

    % Plot spacecraft position at TLI (departure)
    plot3(x(1), y(1), z(1), 'o', ...
        'MarkerEdgeColor', 'k', 'MarkerFaceColor', 'k', 'MarkerSize', 3);

    % Plot spacecraft position at perilune (closest approach)
    plot3(x(imin), y(imin), z(imin), 'o', ...
        'MarkerEdgeColor', 'k', 'MarkerFaceColor', 'k', 'MarkerSize', 2);

    % Plot spacecraft position at end of simulation
    plot3(x(end), y(end), z(end), 'o', ...
        'MarkerEdgeColor', 'r', 'MarkerFaceColor', 'r', 'MarkerSize', 3);

    % Plot Moon at TLI position
    text(xm(1), ym(1), zm(1), 'Moon at TLI');
    Moon1 = surfl(Rm * xx + xm(1), Rm * yy + ym(1), Rm * zz + zm(1));
    set(Moon1, 'FaceAlpha', 0.99);
    shading interp

    % Plot Moon at perilune position
    Moon2 = surfl(Rm * xx + xm(imin), Rm * yy + ym(imin), Rm * zz + zm(imin));
    set(Moon2, 'FaceAlpha', 0.99);
    shading interp

    % Plot Moon at end of simulation
    Moon3 = surfl(Rm * xx + xm(end), Rm * yy + ym(end), Rm * zz + zm(end));
    set(Moon3, 'FaceAlpha', 0.99);
    shading interp

    % Adjust figure visualization settings
    axis image  % Equal scaling for all axes
    axis vis3d  % Fixed aspect ratio for 3D rotation
    axis off    % Remove axis lines and labels
    view([1, 1, 1]);  % Set view angle (isometric)
end
