function [U, V, W, up, vp, wp] = reynoldsDecomposition(u, v, w)
% 
% Matt Werner (m.werner@vt.edu) - June 24, 2021
% 
% Decompose a measured fluid velocity field in the sense of Reynolds
% (average + fluctuating) at a given (fixed) time t. The provided data are
% generally matrix valued, providing description over a volume of positions
% (x, y, z) with increasing index in the matrix corresponding to increasing
% position. The matrix is ordered in the following sense:
%   1. (z, :, :) The elevation  (z) fills the rows
%   2. (:, x, :) The downrange  (x) fills the columns
%   3. (:, :, y) The crossrange (y) fills the "pages".
% (If data is provided in an (x,y,z) format, then `permute(k, [2, 3, 1])'
%  will place them in proper ordering, where k = u, v, w).
% 
% 
% Note: For a wind field assumed to be homogeneous at a given elevation
%       (that is, u = u(z), v = v(z), and w = w(z), i.e. independent of x
%       and y), these quantities (u, v, w) correspond to a instrument tower
%       measuring the wind speed at various (n) elevations. These values
%       are then vector-valued with their FIRST elements corresponding to
%       the measurement made at the LOWEST elevation on the tower. The
%       dimension of the vector is therefore n-by-1(-by-1), the ordinary
%       column vector.
% 
%       Proper ordering is not required for the Reynolds decomposition
%       itself, but is required when gradients need to be evaluated. Thus,
%       it is best to keep the data ordered from the beginning.
% 
%    Inputs:
% 
%           u, v, w - Downwind (u), crosswind (v), and vertical (w) wind
%                     velocities at given locations (x, y, z) and time t
%                     on which Reynolds decomposition is to be performed.
%                     The decomposition provides separation of mean and
%                     fluctuating components as
%                                 k = mean(k) + (k - mean(k)),
%                                     \_____/   \___________/
%                                        K           k'
%                     where k = u, v, w. Here, K is the average (mean)
%                     value of k and k' = (k - K) is the "fluctuation of
%                     k".
%                     Size: z-by-x-by-y (matrix)
%                     Units: m/s (meters per second)
% 
%    Outputs:
% 
%           U, V, W - The respective average (mean) values of the velocity
%                     components u, v, and w.
%                     Size: z-by-x-by-y (matrix)
%                     Units: m/s (meters per second)
% 
%        up, vp, wp - The respective fluctiation values of the velocity
%                     components u, v, and w.
%                     Size: z-by-x-by-y (matrix)
%                     Units: m/s (meters per second)
% 

% Compute the mean value of each component
U = mean(u);
V = mean(v);
W = mean(w);

% Compute the fluctuation of each component
up = u - U;
vp = v - V;
wp = w - W;