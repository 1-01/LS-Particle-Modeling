function L = ObukhovLength(ustar, T0, AVGwpTp)
% 
% Matt Werner (m.werner@vt.edu) - July 1, 2021
% 
% Calculate the Obukhov length (L), determining the stability of the
% atmosphere. This quantity represents "the ratio of mechanically-derived
% total kinetic energy (TKE) to buoyancy-derived TKE" (Aylor, Eq. 8.11).
% 
% Given an elevation z, atmospheric conditions are considered...
%  NEUTRAL   if   z/|L| << 1
% UNSTABLE   if   L < 0  AND  z/|L| > 0.1
%   STABLE   if   L > 0  AND  z/|L| > 0.1.
% 
%    Inputs:
% 
%             ustar - Friction velocity.
%                     Size: 1-by-1 (scalar)
%                     Units: m/s (meters per second)
% 
%                T0 - Potential temperature.
%                     Size: 1-by-1 (scalar)
%                     Units: K (Kelvin)
% 
%           AVGwpTp - Time-averaged value of the vertical wind fluctuation
%                     with the temperature fluctuation. This quantity is
%                     identical to the product of the local sensible heat
%                     flux (H0) with the inverse of air density (rho) and
%                     specific heat at constant pressure (cp), or
%                                   avgwpTp = H0/(rho*cp).
%                     Size: 1-by-1 (scalar)
%                     Units: K m/s (Kelvin * meters per second)
% 
%    Outputs:
% 
%                 L - The Obukhov length.
%                     Size: 1-by-1 (scalar)
%                     Units: m (meters)
% 

%% Checks
% No checks

%% Computation
% Define constants
k = 0.4; % von Karman constant
g = 9.81; % gravitational acceleration (m/s2)
% Calculate the Obukhov length
L = -ustar^3 / (k*(g/T0)*AVGwpTp);