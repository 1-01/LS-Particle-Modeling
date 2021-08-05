function windDirectionTo = convertDirectionTailToHead(windDirectionFrom)
% 
% Matt Werner (m.werner@vt.edu) - August 4, 2021
% 
% Calculate the angle in which the wind is blowing, rather than the angle
% where the wind is coming from. 
% 
% That is, given the angle from where the wind is coming from (x), the
% angle describing the direction of where the wind is going (y) is
%                       y = (x + 180) mod 360.
% Here, x is measured FROM due North in the clockwise sense such that...
% 
%   1. x = 0   means wind is blowing FROM the north
%      y = 180 means wind is blowing TOWARDS the south
% 
%   2. x = 90  means wind is blowing FROM the east
%      y = 270 means wind is blowing TOWARDS the west
% 
%   3. etc.
% 
%    Inputs:
% 
% windDirectionFrom - The direction of the wind, measured clockwise from
%                     due North, describing where the wind is coming FROM.
%                     Size: n-by-1 (vector)
%                     Units: deg (degrees)
% 
%   windDirectionTo - The direction of the wind, measured clockwise from
%                     due North, describing where the wind is going TOWARDS.
%                     Size: n-by-1 (vector)
%                     Units: deg (degrees)
% 

%% Checks
% No checks

%% Computation
% Swing the direction of the wind around by 180 degrees to describe where
% the wind is blowing towards rather than coming from.
windDirectionTo = windDirectionFrom + 180;

% Modulate the wind direction to be between 0 and 360 degrees
windDirectionTo = mod(windDirectionTo, 360);