function [jd] = julian_day(year, month, day, UT)
%==========================================================================
% julian_day: Computes the Julian Day Number (JD) for a given calendar date 
%             and Universal Time (UT).
%
% Inputs:
%   year   - Year (e.g., 2025)
%   month  - Month (1 to 12)
%   day    - Day of the month
%   UT     - Universal Time in hours (decimal format)
%
% Output:
%   jd     - Julian Day Number (JD), including the fractional day.
%
% Formula Reference:
%   Based on the standard astronomical algorithm for Julian Date conversion.
%==========================================================================

    % Compute the Julian Day at 0h UT (midnight)
    j0 = 367 * year ...                                  % Year term
        - fix(7 * (year + fix((month + 9) / 12)) / 4) ...% Century correction term
        + fix(275 * month / 9) ...                       % Month term
        + day ...                                        % Day term
        + 1721013.5;                                     % Julian Day offset to start at J2000 epoch

    % Add the fractional day corresponding to Universal Time (UT in hours)
    jd = j0 + UT / 24;

end
