function [r0, v0, y0, w0_hat] = launchConditions2RV(Re, dec0, alt0, RAAN0, gamma0, v0_mag, rm0)
%==========================================================================
% launchConditions2RV: Computes the initial position and velocity vectors
%                      (r0, v0) of the spacecraft at Trans-Lunar Injection (TLI),
%                      given the departure conditions.
%
% Inputs:
%   Re      - Radius of the Earth [km]
%   dec0    - Declination angle of departure point [degrees]
%   alt0    - Altitude above Earth's surface [km]
%   RAAN0   - Right Ascension of Ascending Node (RAAN) [degrees]
%   gamma0  - Flight Path Angle at TLI [degrees]
%   v0_mag  - Magnitude of velocity at TLI [km/s]
%   rm0     - Moon's position vector at TLI in ECI frame [km]
%
% Outputs:
%   r0      - Initial position vector of the spacecraft in ECI [km]
%   v0      - Initial velocity vector of the spacecraft in ECI [km/s]
%   y0      - Initial state vector [r0_x, r0_y, r0_z, v0_x, v0_y, v0_z]
%   w0_hat  - Unit vector normal to the plane formed by r0 and rm0 (binormal)
%==========================================================================

    % Convert angles from degrees to radians for computation
    dec_rad = deg2rad(dec0);     % Declination [rad]
    RAAN_rad = deg2rad(RAAN0);   % RAAN [rad]
    gamma_rad = deg2rad(gamma0); % Flight Path Angle [rad]

    % Compute the orbital radius at departure (Earth's radius + altitude)
    r0_mag = alt0 + Re;  % [km]

    % Compute components of initial position vector (r0) in ECI frame
    r0_x = r0_mag * cos(dec_rad) * cos(RAAN_rad);  % X component [km]
    r0_y = r0_mag * cos(dec_rad) * sin(RAAN_rad);  % Y component [km]
    r0_z = r0_mag * sin(dec_rad);                  % Z component [km]

    % Assemble position vector r0
    r0 = [r0_x, r0_y, r0_z];  % [km]

    % Compute unit radial vector (ur_hat)
    ur_hat = r0 / norm(r0);

    % Compute unit normal vector to the orbital plane (w0_hat)
    w0_hat = cross(r0, rm0) / norm(cross(r0, rm0));

    % Compute unit transverse vector (ut_hat), perpendicular to r0 and w0_hat
    ut_hat = cross(w0_hat, ur_hat);

    % Compute initial velocity vector (v0) by combining radial and transverse components
    v0 = v0_mag * sin(gamma_rad) * ur_hat + v0_mag * cos(gamma_rad) * ut_hat;
    % Assemble initial state vector y0 = [r0_x, r0_y, r0_z, v0_x, v0_y, v0_z]
    y0 = [r0_x, r0_y, r0_z, v0];

    % Diagnostic print statements
    % fprintf('Initial position vector of Probe (r0) = [%.2f, %.2f,%.2f] km\n', r0);  % Position vector (r0)
    % fprintf('Initial velocity vector of Probe (v0) = [%.4f, %.4f, %.4f] km/s\n', v0);  % Velocity vector (v0)
    % fprintf('Initial state vector of Probe in ECI frame (y0_prb) = [%.4f, %.4f, %.4f, %.4f, %.4f, %.4f]\n', y0); % State vector
end