function u10 = estimateWindSpeedAt10m(u, z)
% 
% Matt Werner (m.werner@vt.edu) - August 5, 2021
% 
% Estimate the wind speed 10 meters above the ground's surface given a
% measured wind speed from a known height.
% 
% This estimation is done using only one of several (5) different 
% methods. These methods come from (*) and are labeled as follows.
%   1. The Logarithmic Law
%   2. The Power Law (1)       <--- This is the current method being used
%   3. The Power Law (2)
%   4. The Exponential Law (1)
%   5. The Exponential Law (2)
% 
% NOTE: The source (*) discusses estimating wind speed u from 0 to 10m
%       given a wind speed measurement u10 at 10m altitude above the
%       ground. This application, however, simply switches the models
%       around in order to estimate the wind speed at 10 meters given a
%       measurement from 0 to 10m.
% 
% (*) Chen, Bundy, Hoff. Modeling the Variation of Wind Speed with Height
%       for Agricultral Source Pollution Control. 1998.
% 
%    Inputs:
% 
%                 u - The measured wind speed at the known altitude z.
%                     Size: 1-by-1 (scalar)
%                     Units: m/s (meters per second)
% 
%                 z - The altitude at which the given wind speed u is
%                     measured.
%                     Size: 1-by-1 (scalar)
%                     Units: m (meters)
% 
%    Outputs:
% 
%               u10 - The estimated wind speed at an altitude of 10 meters.
%                     Size: 1-by-1 (scalar)
%                     Units: m/s (meters per second)

%% Checks
% Ensure that the measured wind speed was taken at an altitude of less than
% 10 meters.
assert(0 < z & z < 10, "Measured wind speed must be between 0 and 10 meters.")

%% Computation
% Compute the wind speed at 10 meters altitude using the power law (1)
% method. 
% 
% Use the Sutton (1953) estimate of the exponent (1/7)
b = 1/7;
u10 = u*(Z/10)^-b;