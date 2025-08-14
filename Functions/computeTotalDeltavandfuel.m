function [deltaV, fuel_used] = computeTotalDeltavandfuel(alt0, mu_e, fac, vesc, Re, Isp, m0, g0)
%==========================================================================
% computeTotalDeltavandfuel: Calculates the total ΔV required for a 
% Trans-Lunar Injection (TLI) maneuver from a given Low Earth Orbit (LEO)
% altitude, and computes the corresponding fuel consumption using the 
% rocket equation.
%
% Inputs:
%   alt0   - Initial altitude of the spacecraft above Earth (km)
%   mu_e   - Standard gravitational parameter of Earth (km^3/s^2)
%   fac    - Fraction of escape velocity to be used for TLI maneuver
%   vesc   - Local escape velocity at initial orbit altitude (km/s)
%   Re     - Radius of Earth (km)
%   Isp    - Specific impulse of the propulsion system (s)
%   m0     - Initial mass of the spacecraft (kg)
%   g0     - Standard gravity (m/s^2)
%
% Outputs:
%   deltaV     - Total ΔV required for TLI (km/s)
%   fuel_used  - Fuel mass used for the maneuver (kg)
%==========================================================================

    % Compute the orbital radius from Earth's center to the spacecraft (km)
    r0 = Re + alt0;
    
    % Compute the circular orbital velocity at radius r0 (km/s)
    v_circular = sqrt(mu_e / r0);
    
    % Compute the velocity required for TLI as a fraction of the escape velocity (km/s)
    v_TLI = fac * vesc;
    
    % Compute the required delta-V (difference between TLI velocity and circular orbit velocity)
    deltaV = v_TLI - v_circular;
    
    % Compute the fuel mass used using the rocket equation
    fuel_used = computeFuelUsed(deltaV, Isp, g0, m0);
end