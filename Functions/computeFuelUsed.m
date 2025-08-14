function fuel = computeFuelUsed(deltaV, Isp, g0, m0)
%==========================================================================
% computeFuelUsed: Calculates the fuel mass used for a given deltaV using 
%                  the Tsiolkovsky rocket equation.
%
% Inputs:
%   deltaV - Delta-V required for the maneuver [km/s]
%   Isp    - Specific Impulse of the propulsion system [s]
%   g0     - Standard gravitational acceleration [m/s^2]
%   m0     - Initial mass of the spacecraft before the burn [kg]
%
% Output:
%   fuel   - Fuel mass consumed to achieve the deltaV [kg]
%
% If deltaV is zero or negative, the function returns fuel = 0.
%==========================================================================
    
    % Check if deltaV is zero or negative (no burn needed)
    if deltaV <= 0
        fuel = 0;  % No fuel is consumed
        return;    % Exit the function
    end

    % Convert deltaV from km/s to m/s for consistency with Isp and g0
    deltaV_mps = deltaV * 1000; 

    % Compute final mass after burn using Tsiolkovsky rocket equation
    % mf = m0 / exp(deltaV / (Isp * g0))
    mf = m0 / exp(deltaV_mps / (Isp * g0));

    % Compute fuel consumed as the difference between initial and final mass
    fuel = m0 - mf;
end