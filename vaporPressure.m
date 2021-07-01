function [e, es, varargout] = vaporPressure(T, RH)
% 
% Matt Werner (m.werner@vt.edu) - July 1, 2021
% 
% Calculate the actual vapor pressure (e) and saturated vapor pressure (es)
% via the Magnus formula.
%                                   /   A t   \
%                      e  = C * exp|  -------  |
%                       s           \  B + t  /
% The coefficients A, B, C are recommended by Alduchov and Eskridge (1996),
% resulting in a relative error of < 0.4% over the range of temperatures
% -40C < T < 50C)*.
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
%             e, es - Actual/saturated vapor pressure, respectively.
%                     Size: 1-by-1 (scalar)
%                     Units: Pa (Pascals)
% 
%         varargout - Coefficients (A1, B1, C1) determined by Alduchov and
%                     Eskridge (1996) for the Magnus formula.
%                     Size: 1-by-3 (cell)
%                     Units: {-, C, Pa} ({unitless, Celsius, Pascals})

%% Checks
% Check that the temperature is within the range for which the Magnus
% formula is valid.
assert(-40 < T & T < 50, "Temperature (%2.2f C) is out of Magnus range", T)
% Check that the relative humidity is between 0 and 100%
assert(0 < RH & RH < 100, "Relative humidity (%3.1f%) is out of bounds")

%% Computation
% Calculate the saturation vapor pressure according to the Magnus formula
% and coefficients recommended by Alduchov, Eskridge (1996)
A1 = 17.625; % No units
B1 = 243.04; % Celsius
C1 = 610.94; % Pascals
es = C1*exp(A1*T / (B1 + T));

% Calculate the actual water vapor pressure according to the definition of
% relative humidity
e = (RH/100)*es;

% Assign varargout as the coefficients of the model
varargout{1} = A1;
varargout{2} = B1;
varargout{3} = C1;