function wstar = convectiveVelocityScale(BLmixingDepth, AVGvirtualT, Hsurface)
% 
% Matt Werner (m.werner@vt.edu) - July 1, 2021
% 
% Calculate the convective velocity scale w*. This quantity, along with the
% boundary layer mixing depth (z_i in Aylor's literature) serve "as the main
% input parameters for convective boundary layer (CBL) transport models"
% (Aylor 9.4.1).
% 
% The convective velocity scale w* is defined
% 
%                       /  g * BLmixingDepth              \ 1/3
%                 w* = |  -------------------- * Hsurface  |   ,
%                       \     AVGvirtualT                 /
% 
% where g = 9.81 (m/s2) is the typical standard value for gravitational
% acceleration.
% 
% The surface heat flux, Hsurface, is identically equivalent to the
% averaged value of the product of vertical wind fluctuation and virtual
% temperature, both evaluated at the surface of Earth, or

% 
%    Inputs:
% 
%     BLmixingDepth - The boundary layer mixing depth (also called z_i).
%                     Size: 1-by-1 (scalar)
%                     Units: m (meters)
% 
%       AVGvirtualT - The time-averaged virtual temperature.
%                     Size: 1-by-1 (scalar)
%                     Units: K (Kelvin)
% 
%          Hsurface - The surface heat flux. This quantity is identically
%                     equivalent to the time-averaged value of the product
%                     of vertical wind fluctuation (w') and virtual
%                     temperature, both evaluated at the earth's surface.
%                     In other words, 
%                                  __________________
%                       Hsurface =  w' * AVGvirtualT   @ earth's surface.
%                     Size: 1-by-1 (scalar)
%                     Units: K m/s (Kelvin * meters per second)
% 

%% Checks
% No checks

%% Computation
% Calculate the convective velocity scale
wstar = (9.81*BLmixingDepth*Hsurface/AVGvirtualT)^(1/3);