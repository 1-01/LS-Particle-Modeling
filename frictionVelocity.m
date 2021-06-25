function ustar = frictionVelocity(meanupwp, meanvpwp)
% 
% Matt Werner (m.werner@vt.edu) - June 25, 2021
% 
% Calculate the friction velocity, u*. The friction velocity is generally
% defined by (Aerial Dispersion of Pollen and Spores, Aylor, Eq. 7.18b)
%                        / _____ 2    _____ 2 \ 1/4
%                  u* = |  u' w'   +  v' w'    |   ,
%                        \                    /
% where appears the mean product of the fluctuations in the horizontal
% directions with the vertical direction.
% 
% Note: This form of u* is used as it provides a more general description
%       for the orientation of the coordinate system. In the case that the
%       coordinate system has its x axis aligned with the mean flow
%       direction, the mean product of (v' w') vanishes.
% 
%    Inputs:
% 
%      mean(u|v)pwp - Mean product of the fluctuations in u and v
%                     corresponding to the x and y directions,
%                     respectively.
%                     Size: 1-by-1 (scalar)
%                     Units: m2/s2 (square meter per square second)
% 
%    Outputs:
% 
%             ustar - The friction velocity.
%                     Size: 1-by-1 (scalar)
%                     Units: m/s (meters per second)
% 

%% Checks
% No checks

%% Computation
% Calculate the friction velocity
if (nargin == 1)
    % Calculate the friction velocity in the sense that the x axis is
    % aligned exactly with the wind such that there is absolutely no
    % correlation with the crossrange and vertical wind fluctuations
    ustar = sqrt(abs(meanupwp));
else
    % Calculate the friction velocity in the sense that the x axis is not
    % aligned with the wind direction (Stull, 1988)
    ustar = (meanupwp^2 + meanvpwp^2)^0.25;
end