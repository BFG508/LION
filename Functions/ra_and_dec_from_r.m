function [ra, dec] = ra_and_dec_from_r(r)
%==========================================================================
% ra_and_dec_from_r: Computes the Right Ascension (RA) and Declination (Dec)
%                    from a given position vector in the ECI frame.
%
% Inputs:
%   r  - Position vector in ECI coordinates [x; y; z] (km)
%
% Outputs:
%   ra  - Right Ascension in degrees [0, 360)
%   dec - Declination in degrees [-90, 90]
%
% The function computes:
%   - Declination as the arcsine of the z-component over |r|
%   - Right Ascension as the angle in the XY-plane measured from X-axis
%     (corrected depending on the sign of the Y-component)
%==========================================================================

    % Normalize components of position vector
    l = r(1)/norm(r);  % Direction cosine along X-axis
    m = r(2)/norm(r);  % Direction cosine along Y-axis
    n = r(3)/norm(r);  % Direction cosine along Z-axis

    % Compute Declination (Dec) in degrees
    dec = asind(n);

    % Compute Right Ascension (RA) depending on the sign of Y-component
    if m > 0
        ra = acosd(l / cosd(dec));
    else
        ra = 360 - acosd(l / cosd(dec));
    end
end
