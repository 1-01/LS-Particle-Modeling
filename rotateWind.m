function [u, v] = rotateWind(N, E, windDirection)
% 
% Matt Werner (m.werner@vt.edu) - June 29, 2021
% 
% Rotate the components of the velocity field describing the wind by the
% angle WINDDIRECTION. WINDDIRECTION is the angle subtended clockwise from
% due North describing where the wind is coming from.
%   - A wind direction of 00 indicates that the wind is coming from the
%     North (i.e. blowing to the South)
%   - A wind direction of 90 indicates that the wind is coming from the
%     East  (i.e. blowing to the West)
%   - etc.
% Note: The wind direction (WINDDIRECTION) is a scalar for the entire flow
%       field. It corresponds to the average, or mean, flow direction over
%       the entire field. Thus, WINDDIRECTION represents the angle made
%       between the East/North/Vertical frame and the
%       downwind/crosswind/vertical frame.
% 
% Supplying the wind components in the East/North/Vertical reference frame
% and the wind direction as defined above results in components describing
% the wind in the downwind (u) and crosswind (v) directions.
% 
%    Inputs:
% 
%              N, E - The North and East components of the wind at each
%                     location in the vector field. These two elements are
%                     paired together and rotated (pointwise) so that the
%                     resulting vector field points in the direction of the
%                     wind.
%                     Size: n-by-m (matrix)
%                     Units: m/s (meters per second)
% 
%     windDirection - Angle measured clockwise from due North describing
%                     the direction of the main movement of wind
%                     Size: 1-by-1 (scalar)
%                     Units: deg (degrees)
% 
%    Outputs:
% 
%              u, v - The downwind and crosswind comopnents of the wind at
%                     each location in the vector field.
%                     Size: n-by-m (matrix)
%                     Units: m/s (meters per second)
% 

%% Checks
% No checks

%% Computation
% Add an additional 180 degrees rotation onto the wind direction to
% describe where the wind is going rather than where it's coming from
windDirection = windDirection + 180;

% Rotate the vector field pointwise into the mean direction of the wind.
% Because this direction is constant, the rotation can occur by expressing
% the rotation matrix as an algebraic sum acting on each element of the
% field.
u = +N*cosd(windDirection) + E*sind(windDirection);
v = -N*sind(windDirection) + E*cosd(windDirection);