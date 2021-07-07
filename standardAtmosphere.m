function [T, p, d] = standardAtmosphere(z, Tsurface, psurface)
% 
% Matt Werner (m.werner@vt.edu) - July 7, 2021
% 
% Calculate the US 1976 Standard Atmosphere temperature for a given
% elevation (z) above ground level (AGL) beneath the tropopause. This
% implementation of the standard atmosphere is normally valid only to 11 km
% altitude (using geometric altitude above sea level (ASL)), in which case
% the lapse rate is
%                  a = -0.0065 K/m,      (0 < HASL < 11000 m).
% Because the elevation about ground level is used, however, the range over
% which the model is defined is modified such that
%                  a = -0.0065 K/m,      (0 < zAGL < 10000 m).
% 
% Note: This implementation of the standard atmosphere uses the measured
%       surface temperature as the reference and altitude is measured
%       relative to the ground, NOT sea level.
% 
% Note: The standard atmosphere calculates its parameters relative to the
%       geopotential (not geometric) altitude. However, even at 45,000 ft
%       (14 km, above the tropopause), the difference between the two is
%       only about 100 ft. As such, the conversion from geometric to
%       geopotential is ignored and the geometric height is used in its
%       place.
% 
%    Inputs:
% 
%                 z - Geometric altitude above ground level.
%                     Size: 1-by-1 (scalar)
%                     Units: m (meters)
% 
%          Tsurface - Measured temperature at the Earth's surface.
%                     Size: 1-by-1 (scalar)
%                     Units: K (Kelvin)
% 
%          psurface - Measured pressure at the Earth's surface.
%                     Size: 1-by-1 (scalar)
%                     Units: Pa (Pascals)
% 
%    Outputs:
% 
%                 T - "Standard" atmospheric temperature.
%                     Size: 1-by-1 (scalar)
%                     Units: K (Kelvin)
% 
%                 p - "Standard" atmospheric pressure.
%                     Size: 1-by-1 (scalar)
%                     Units: Pa (Pascals)
% 
%                 d - "Standard" atmospheric density.
%                     Size: 1-by-1 (scalar)
%                     Units: kg/m3 (kilograms per cubic meter)
% 

%% Checks
% Ensure that the elevation is within 1 km below ground level and 10 km
% above ground level. Note that this quantity is not the actual elevation
% (above sea level), but rather the elevation above ground level. Note that
% an elevation below ground level could correspond to a cliff drop-off or
% similar scenario.
assert(-1000 < z & z < 10000)

%% Computations
% Define the lapse rate, gravity at sea level, and specific gas constant
% for DRY air (no humidity)
dTdz = -0.0065; % K/m
g0 = 9.80665; % m/s2
R = 287.058; % J/(kg K)

% Calculate air density at the surface using the ideal gas law
dsurface = psurface/(R*Tsurface);

% Calculate the exponent used for evaluating the pressure and density
negative_g0OVERdTdzR = -g0/(dTdz*R);

% Calculate the "standard" temperature, pressure, and density at the
% elevation z
T = Tsurface + dTdz*z;
p = psurface*(T/Tsurface)^negative_g0OVERdTdzR;
d = dsurface*(p/psurface)/(T/Tsurface);