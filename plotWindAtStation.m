function plotWindAtStation(z, u, v, w)
% 
% Matt Werner (m.werner@vt.edu) - June 28, 2021
% 
% Plot the wind flow field measured at the station.
% 
%    Inputs:
% 
%                 z - Elevations at which the wind is measured.
%                     Size: n-by-1
%                     Units: m (meters)
% 
%           u, v, w - Velocity components of the wind. These components may
%                     be either in the East/North/Vertical reference frame
%                     or in the downrange/crossrange/vertical reference
%                     frame.
%                     Size: n-by-1
%                     Units: m/s (meters per second)
% 
%    Outputs:
% 
%                   -
% 

% Plot the velocity field of the wind located at the origin
quiver3(zeros(size(z)), zeros(size(z)), z, u, v, w)
axis equal

% Scale the axis so the x-y plane is square with the origin in the middle
scale = max([xlim, ylim]);
xlim([-scale, scale])
ylim([-scale, scale])