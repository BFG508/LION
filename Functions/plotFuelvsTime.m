function plotFuelvsTime(t, m0, deltaV_total, burn_duration, Isp)
%==========================================================================
% plotFuelvsTime: Plots the cumulative fuel consumption profile during 
%                 the main maneuver burn, assuming a constant rate burn.
%
% Inputs:
%   t             - Time vector from ODE integration [s]
%   m0            - Initial mass of the spacecraft [kg]
%   deltaV_total  - Total ΔV applied during the burn [km/s]
%   burn_duration - Duration of the main burn [s]
%   Isp           - Specific impulse of the engine [s]
%
% The function computes the fuel used over time, with a linear increase
% during the burn and constant afterwards, and plots the result.
%==========================================================================

    g0 = 9.80665; % Standard gravity [m/s^2]

    % Convert deltaV from km/s to m/s
    deltaV_mps = deltaV_total * 1000;  % [m/s]

    % Compute total fuel used with the rocket equation
    fuel_total = m0 * (1 - exp(-deltaV_mps / (Isp * g0)));  % [kg]

    % Initialize fuel profile array with zeros (same size as time vector t)
    fuel_profile = zeros(size(t));

    % Logical index vector: true for times during burn, false afterwards
    idx_burn = t <= burn_duration;

    % Linearly distribute fuel usage during the burn duration
    fuel_profile(idx_burn) = linspace(0, fuel_total, sum(idx_burn));

    % After burn completion, fuel remains constant at total fuel used
    fuel_profile(~idx_burn) = fuel_total;

    % Plot the fuel consumption profile vs time (in hours)
    figure;
    plot(t/3600, fuel_profile, 'r', 'LineWidth', 2);
    xlabel('Time (h)', 'Interpreter', 'latex');
    ylabel('Fuel used (kg)', 'Interpreter', 'latex');
    title('\bf{Main maneuver fuel usage}', 'Interpreter', 'latex');
    grid on;
end