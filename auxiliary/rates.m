function dydt = rates(t, y, jd0, ttt, days, mu_m, mu_e)
%==========================================================================
% rates: Computes the derivatives of the state vector (position and velocity)
%        for a spacecraft influenced by Earth's and Moon's gravity.
%        This function is designed to be used with an ODE solver (e.g., ode45).
%
% Inputs:
%   t    - Current time (s)
%   y    - Current state vector [X; Y; Z; vX; vY; vZ]
%   jd0  - Julian Date of TLI (Trans-Lunar Injection)
%   ttt  - Target time to reach the Moon (s)
%   days - Conversion factor (1 day in seconds)
%   mu_m - Gravitational parameter of the Moon (km^3/s^2)
%   mu_e - Gravitational parameter of the Earth (km^3/s^2)
%
% Outputs:
%   dydt - Derivative of the state vector [vX; vY; vZ; aX; aY; aZ]
%==========================================================================

    % Compute Julian Date at current time t
    jd = jd0 - (ttt - t) / days;

    % Extract position and velocity components from state vector y
    X = y(1);
    Y = y(2);
    Z = y(3);
    vX = y(4);
    vY = y(5);
    vZ = y(6);

    % Compute spacecraft position vector and its magnitude
    r = [X; Y; Z];
    r_mag = norm(r);

    % Compute Moon's position at current Julian Date using ephemeris
    [rm, ~] = simpsons_lunar_ephemeris(jd);
    rm_mag = norm(rm);

    % Vector from spacecraft to Moon (r_moon - r_spacecraft)
    rms = rm - r;
    rms_mag = norm(rms);

    % Compute acceleration due to Earth's gravity (central body term)
    a_earth = -mu_e * r / r_mag^3;

    % Compute acceleration due to Moon's gravity (including indirect term)
    a_moon = mu_m * (rms / rms_mag^3 - rm / rm_mag^3);

    % Total acceleration acting on the spacecraft
    a = a_earth + a_moon;

    % Extract acceleration components
    aX = a(1);
    aY = a(2);
    aZ = a(3);

    % Assemble derivative of the state vector (velocity and acceleration)
    dydt = [vX; vY; vZ; aX; aY; aZ];

    % Diagnostic print statements
    % fprintf('y = [%.4f, %.4f, %.4f, %.4f, %.4f, %.4f]\n', y);
    % fprintf('dydt = [%.4f, %.4f, %.4f, %.4f, %.4f, %.4f]\n', dydt);
end
