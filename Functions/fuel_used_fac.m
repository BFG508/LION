function fuel = fuel_used_fac(fac, alt0, mu_e, vesc, Re, Isp, m0, g0)
%==========================================================================
% fuel_used_fac: Computes the amount of fuel used for a given scaling factor 
%                'fac' applied to the escape velocity.
%
% Inputs:
%   fac    - Scaling factor applied to escape velocity (vesc)
%   alt0   - Initial altitude above Earth's surface [km]
%   mu_e   - Gravitational parameter of Earth [km^3/s^2]
%   vesc   - Escape velocity at the initial orbit [km/s]
%   Re     - Earth's radius [km]
%   Isp    - Specific Impulse of the propulsion system [s]
%   m0     - Initial mass of the spacecraft before burn [kg]
%   g0     - Standard gravitational acceleration [m/s^2]
%
% Output:
%   fuel   - Fuel mass consumed to achieve the required deltaV [kg]
%
% This function is typically used as an objective function in optimization.
%==========================================================================

    % Compute orbital radius (distance from Earth's center)
    r0 = Re + alt0;  % [km]

    % Compute circular orbit velocity at altitude alt0 (LEO speed)
    v_circular = sqrt(mu_e / r0);  % [km/s]

    % Compute Translunar Injection (TLI) velocity by scaling escape velocity
    v_TLI = fac * vesc;  % [km/s]

    % Compute required deltaV for the maneuver
    deltaV = v_TLI - v_circular;  % [km/s]

    % Compute fuel consumed using the deltaV and rocket equation
    fuel = compute_fuel_used(deltaV, Isp, g0, m0);  % [kg]
end
