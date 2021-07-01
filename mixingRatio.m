function [w, ws] = mixingRatio(T, RH, p)
% 
% Matt Werner (m.werner@vt.edu) - July 1, 2021
% 
% Calculate the actual/saturated mixing ratios. The actual mixing ratio (w)
% is given by
%                                         e
%                           w = 621.97 -------,
%                                       p - e
% whereas the saturated mixing ratio (ws) is defined in terms of the
% relative humidity (RH, in percent) as
% 
%                           RH = 100 * w/ws.
% 
%    Inputs:
% 
%                 T - Ambient temperature.
%                     Size: 1-by-1 (scalar)
%                     Units: C (degrees Celsius)
% 
%                RH - Relative humidity.
%                     Size: 1-by-1 (scalar)
%                     Units: % (percent)
% 
%                 p - Ambient pressure.
%                     Size: 1-by-1 (scalar)
%                     Units: Pa (Pascals)
% 
%    Outputs:
% 
%             w, ws - Actual and saturated mixing ratios, respectively.
%                     Size: 1-by-1 (scalar)
%                     Units: g/g (grams per gram)
% 

%% Checks
% No checks

%% Computation
% Obtain the actual vapor pressure
e = vaporPressure(T, RH);

% Calculate the mixing ratio
w = 621.97 * e/(p - e); % (kg/g)
% Calculate the saturated mixing ratio by the definition of relative
% humidity
ws = w/(RH/100); % (kg/g)

% Convert the mixing ratios to (g/g)
w = w/100; ws = ws/100;