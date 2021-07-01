function virtualT = virtualTemperature(T, RH, p, p0)
% 
% Matt Werner (m.werner@vt.edu) - July 1, 2021
% 
% Calculate the virtual temperature according to the standard atmosphere.
% The virtual temperature is, in general, different from the potential
% (total) temperature as it accounts for humidity in the air. The
% correction is given by
% 
%           virtualT  = potentialT * (1 + 0.61 r),
% 
% where r is the mixing ratio of the air.
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
%                p0 - Static pressure (pressure measured when the wind is
%                     not blowing).
%                     Size: 1-by-1 (scalar)
%                     Units: Pa (Pascals)
% 
%    Outputs:
% 
%          virtualT - Virtual temperature.
%                     Size: 1-by-1 (scalar)
%                     Units: C (degrees Celsius)
% 

%% Checks
% No checks

%% Computation
% Calculate the potential temperature using the standard atmosphere for
% elevations within the troposphere
R_over_cp_dryAir = 0.286;
T0 = T*(p0/p)^R_over_cp_dryAir;

% Correct the potential temperature to account for relative humidity (i.e.
% having moisture in the air, air that is not completely dry)
r = mixingRatio(T, RH, p);
virtualT = T0*(1 + 0.61*r);