function [pos, vel] = simpsonsLunarEphemeris(jd)
%==========================================================================
% simpsonsLunarEphemeris: Computes an approximate position and velocity
% vector of the Moon at a given Julian Date using a simplified ephemeris
% based on harmonic terms.
%
% Inputs:
%   jd  - Julian Date (days)
%
% Outputs:
%   pos - Position vector of the Moon in ECI frame [x; y; z] (km)
%   vel - Velocity vector of the Moon in ECI frame [vx; vy; vz] (km/s)
%==========================================================================

    % Conversion factor: seconds in a Julian century
    tfac = 36525 * 3600 * 24;  % s/century

    % Time in Julian centuries since J2000.0 epoch (Jan 1, 2000, 12:00 TT)
    t = (jd - 2451545.0) / 36525;

    % Amplitude coefficients (in km) for harmonic series (3 components x 7 terms)
    a = [383.0 31.5 10.6 6.2 3.2 2.3 0.8;
         351.0 28.9 13.7 9.7 5.7 2.9 2.1;
         153.2 31.5 12.5 4.2 2.5 3.0 1.8] * 1e3;

    % Frequencies (in radians per century) for each term
    b = [8399.685 70.990 16728.377 1185.622 7143.070 15613.745 8467.263;
         8399.687 70.997 8433.466 16728.380 1185.667 7143.058 15613.755;
         8399.672 8433.464 70.996 16728.364 1185.645 104.881 8399.116];

    % Phase angles (in radians) for each term
    c = [5.381 6.169 1.453 0.481 5.017 0.857 1.010;
         3.811 4.596 4.766 6.165 5.164 0.300 5.565;
         3.807 1.629 4.595 6.162 5.167 2.555 6.248];

    % Initialize position and velocity vectors
    pos = zeros(3, 1);  % Position vector [x; y; z] (km)
    vel = zeros(3, 1);  % Velocity vector [vx; vy; vz] (km/s)

    % Compute position and velocity contributions from harmonic series
    for i = 1:3  % Loop over x, y, z components
        for j = 1:7  % Loop over each term in the series
            angle = b(i, j) * t + c(i, j);  % Compute phase angle for term (i,j)
            
            % Add sine contribution to position (in km)
            pos(i) = pos(i) + a(i, j) * sin(angle);

            % Add cosine contribution to velocity (proportional to frequency)
            vel(i) = vel(i) + a(i, j) * cos(angle) * b(i, j);
        end
        % Convert velocity from km/century to km/s
        vel(i) = vel(i) / tfac;
    end
end