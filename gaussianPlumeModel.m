function [c, maxcpos] = gaussianPlumeModel(h, u, x, y, z, stability)
% 
% Matt Werner (m.werner@vt.edu) - July 30, 2021
% 
% Calculate the concentration of a Gaussian plume dispersion model. 
% 
% The concentration from this Gaussian model is able to provide an estimate
% as to where the best positioning for a measurement in an experiment might
% be. (Here, 'best positioning' refers to the location providing the
% highest probability of capturing particles.)
% 
% Assumptions of the Gaussian plume are:
%   1. The plume starts from a mathematical point (a point source)
%   2. The source of pollution (Q) is constant
%   3. Wind direction and wind speed are constant in space and time
%   4. Atmospheric turbulence is constant in space and time
%   5. The plume behaves as if it reflects off of the ground
% (Visscher, Air Dispersion Modeling Sec. 2.3.1)
% 
% Because the resulting concentration of the Gaussian plume is directly
% proportional to the release rate Q, we take Q = 1 without loss of
% generality for providing the location (NOT the value) of highest
% probability of capturing particles.
% 
% For stability classes, the terrain is assumed to be RURAL.
% 
%    Inputs:
% 
%                 h - Effective source height.
%                     Size: 1-by-1 (scalar)
%                     Units: m (meters)
% 
%                 u - Wind speed at the effective source height.
%                     Size: 1-by-1 (scalar)
%                     Units: m/s (meters per second)
% 
%           x, y, z - Locations in the downwind, cross-wind and vertical
%                     directions, respectively, at which to evaluate the
%                     Gaussian plume concentration.
%                     Size: n-by-m-by-k (3D array)
%                     Units: m (meters)
% 
%         stability - Stability class of the weather for determining the
%                     dispersion parameters according to the Briggs
%                     equations (1973). These classes are:
%                           A - Very unstable
%                           B - Moderately unstable
%                           C - Slightly unstable
%                           D - Neutral
%                           E - Slightly stable
%                           F - Stable
%                     Note: These stability classes are cAsE-SeNsItIvE.
%                     Note: More sophisticated approaches for determining
%                           these dispersion parameters are available; see
%                           Visccher, Air Dispersion Modeling Sec 6.3 & 6.6.
%                     Size: 1-by-1 (char)
%                     Units: N/A
% 
%    Output:
% 
%                 c - The concentration from the Gaussian plume.
%                     Size: n-by-m-by-k (3D array)
%                     Units: g/m3 (grams per cubic meter)
% 
%           maxcpos - The location of the greatest concentration according
%                     to the given resolution as an (x, y, z) tuple.
%                     Size: 1-by-3 (array)
%                     Units: m (meters)
% 

%% Checks
% Ensure that the effective source height and wind speed are both positive
assert(h > 0, "Effective source height (%2.2f) must be positive scalar.", h)
assert(u > 0, "Wind speed (%2.2f) must be positive scalar.", u)
% Ensure that the stability class is valid (A - F)
assert(length(stability) == 1 && contains('ABCDEF', stability), ...
    "Please provide a stability class (A, B, C, D, E, or F) depending on the weather conditions.")

%% Computation
% Determine the dispersion parameters
switch stability
    case 'A'
        sigmay = 0.22*x.*(1 + 0.0001*x).^-0.5;
        sigmaz = 0.2*x;
    case 'B'
        sigmay = 0.16*x.*(1 + 0.0001*x).^-0.5;
        sigmaz = 0.12*x;
    case 'C'
        sigmay = 0.11*x.*(1 + 0.0001*x).^-0.5;
        sigmaz = 0.08*x.*(1 + 0.0002*x).^-0.5;
    case 'D'
        sigmay = 0.08*x.*(1 + 0.0001*x).^-0.5;
        sigmaz = 0.06*x.*(1 + 0.0015*x).^-0.5;
    case 'E'
        sigmay = 0.06*x.*(1 + 0.0001*x).^-0.5;
        sigmaz = 0.03*x.*(1 + 0.0003*x).^-1;
    case 'F'
        sigmay = 0.04*x.*(1 + 0.0001*x).^-0.5;
        sigmaz = 0.016*x.*(1 + 0.0003*x).^-1;
end

% Calculate the Gaussian plume concentration using the above assumptions
% and equation (Visscher, Eq. 2.2). Note that the emission/release rate Q
% is set to 1, as the location, not the value, of highest concentration is
% what's important.
%
% Also perform the calculation in an efficient manner. This is done by
% utilizing (overriding) the variables already defined for the dispersion
% parameters as the size of the objects is not expected to change.
% > This is being done to minimize computation time in case a massive grid
%   of (x, y, z) is requested.
% > The resulting formula will appear to be incorrect, but it is, in fact,
%   identical to the one provided by (2.2).
c = 1/(2*pi*u*sigmay.*sigmaz);
% Override the previously defined dispersion parameter (for z) since the
% remainder of the formula uses it twice and only in the denominator
sigmay = exp(-0.5*(y./sigmay).^2);
sigmaz = -0.5*sigmaz.^-2;
c = c.*sigmay.*(exp(sigmaz.*(z - h).^2) + exp(sigmaz.*(z + h).^2));

% Determine the location of maximum concentration using its linear index
[~, lind] = max(c, [], 'all', 'linear');
maxcpos = [x(lind), y(lind), z(lind)];