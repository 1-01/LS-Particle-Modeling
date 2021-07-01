function potentialT = potentialTemperature(T, p)
% 
% Matt Werner (m.werner@vt.edu) - July 1, 2021
% 
% Calculate the potential (total) temperature using the standard atmosphere
% model (dry air) within the troposphere. The potential temperature in this
% case is expressed
% 
%                           /  P   /  \  R/cp
%                   T  = T |    0 /    |     ,
%                    0      \    / P  /
% 
% where R/cp ~ 0.286 for dry air. The quantity P0 is the potential (total)
% pressure corresponding to the potential temperature and is taken to be
% its reference value at sea level of exactly 101325 Pa.
% 
%    Inputs:
% 
%                 T - Ambient temperature.
%                     Size: 1-by-1 (scalar)
%                     Units: C (degrees Celsius)
% 
%                 p - Ambient pressure.
%                     Size: 1-by-1 (scalar)
%                     Units: Pa (Pascals)
% 
%    Outputs:
% 
%        potentialT - Potential (total) temperature according to the
%                     standard atmosphere model evaluated within the
%                     troposphere. This quantity corresponds to the
%                     measured values of the ambient temperature (T) and
%                     pressure (p).
%                     Size: 1-by-1 (scalar)
%                     Units: C (degrees Celsius)
% 

% Calculate the potential temperature
potentialT = T*(101325/p)^0.286;