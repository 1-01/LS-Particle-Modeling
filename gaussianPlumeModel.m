function [c, maxcpos, maxcpos1m] = gaussianPlumeModel(h, u, x, y, z, stability)
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
%                          ([x, y, z] = ndgrid(x, y, z)
%                                        OR
%                           1-by-n (x is a row vector)
%                           1-by-1 (y is scalar - probably y = 0)
%                           m-by-1 (z is a column vector)
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
%         maxcpos1m - The location of the greatest concentration nearest to
%                     a 1 meter altitude according to the given resolution
%                     as an (x, y, z) tuple.
%                     Size: 1-by-3 (array)
%                     Units: m (meters)
% 

%% Checks
% Ensure that the effective source height and wind speed are both positive
assert(h > 0, "Effective source height (%2.2f) must be positive scalar.", h)
assert(u > 0, "Wind speed (%2.2f) must be positive scalar.", u)
% Check that (x,y,z) dimensions are appropriate. No issues if they are all
% the same size. Otherwise, x should be a row vector, y a scalar, and z a
% column vector
arrays3D = nan;
if (all([size(x) == size(y), size(y) == size(z)]))
    arrays3D = true;
elseif ((numel(size(x)) == 2 && size(x,1) == 1 && ...
        numel(size(y)) == 2 && isscalar(y) && ...
        numel(size(z)) == 2 && size(z,2) == 1))
    arrays3D = false;
end
assert(~isnan(arrays3D), "Invalid position dimensions (if z is 1-by-m, try z').")
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
c = (1./(2*pi*u*sigmay.*sigmaz)).*exp(-0.5*(y./sigmay).^2)...
    .*(exp(-0.5*((z - h)./sigmaz).^2) + exp(-0.5*((z + h)./sigmaz).^2));

% Find the positions of maximum concentration throughout the field and 1
% meter above the ground. This procedure can differ if the inputs are 3D
% arrays or not.
if (arrays3D)
    % All arrays are the same size. Thus, a linear index will work.
    % Determine the position of maximum concentration using a linear index
    [~, lind] = max(c, [], 'all', 'linear');
    maxcpos = [x(lind), y(lind), z(lind)];
    % Determine the index of maximum concentration at a 1 meter elevation.
    % Since z comes from ndgrid, only check the first element (1,1,:)
    % because all elements in a given stack are equal, i.e. for any i,j in
    % the stack k, (i,j,k) = (:,:,k).
    [~, pos1m] = min(abs(z(1,1,:) - 1));
    [~, lind1m] = max(c(:,:,pos1m), [], 'all', 'linear');
    [row, col] = ind2sub(size(c(:,:,pos1m)), lind1m);
    maxcpos1m = [x(row, col, pos1m), y(row, col, pos1m), z(row, col, pos1m)];
else
    % The concentration is an m-by-n matrix whose maximum value location
    % at (k,i) corresponds to the positions z(k) and x(i), evaluated at the
    % given scalar y.
    [~, I] = max(c, [], 'all', 'linear');
    [row, col] = ind2sub(size(c), I);
    maxcpos = [x(col), y, z(row)];
    % Find the index of z that is closest to 1 meter above the ground and
    % then find the maximum concentration at this altitude
    [~, pos1m] = min(abs(z - 1));
    [~, ind1m] = max(c(pos1m, :));
    maxcpos1m = [x(ind1m), y, z(pos1m)];
end