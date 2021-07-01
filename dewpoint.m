function dewpointT = dewpoint(T, RH)
% 
% Matt Werner (m.werner@vt.edu) - July 1, 2021
% 
% Calculate the dewpoint temperature given the measured ambient temperature
% (T, in Celsius) and relative humidity (RH, in percent).
% 
% The dewpoint temperature is calculated using the Magnus formula, with
% coefficients recommended by Alduchov and Eskridge (1996), to evaluate the
% saturation vapor pressure (which has a relative error of < 0.4% over the
% range of temperatures -40C < T < 50C)*.
% 
% * [Lawrence, M. "The Relationship between Relative Humidity and the
%       Dewpoint Temperature in Moist Air". 
%       https://journals.ametsoc.org/view/journals/bams/86/2/ ...
%           ... bams-86-2-225.xml?tab_body=pdf]
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
%    Outputs:
% 
%         dewpointT - Dewpoint temperature.
%                     Size: 1-by-1 (scalar)
%                     Units: C (degrees Celsius)
% 

%% Checks
% No checks

%% Computation
% Calculate the actual vapor pressure
[e, ~, A1, B1, C1] = vaporPressure(T, RH);
% Calculate the dewpoint temperature
dewpointT = B1*log(e/C1) / (A1 - log(e/C1));