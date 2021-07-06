function epsilon = dTKEdt(z, ustar, wstar, BLmixingDepth, L)
% 
% Matt Werner (m.werner@vt.edu) - July 1, 2021
% 
% Calculate the dissipation rate of the wind's turbulent kinetic energy
% (TKE) within the atmosphere's convective boundary layer (CBL). This 
% expression provides a profile "for which both buoyancy- and
% shear-generated turbulence are important" (Aylor, eq. 9.45).
%                      3                                        -1/4
%                0.4 w*        3  /       z    \    (1 - 15 z/L)
%     epsilon = ---------- + u*  |  1 - ------  | -----------------  ,
%                  z_i            \      z_i   /         k z
% 
% where wstar is the convective velocity scale, ustar is the friction
% velocity, z_i is the boundary layer mixing depth, z is the altitude, and
% L is the Obukhov length. The parameter k here (= 0.4) is von Karman's
% constant.
% 
%    Inputs:
% 
%                 z - Altitude.
%                     Size: 1-by-1 (scalar)
%                     Units: m (meters)
% 
%             ustar - Friction velocity.
%                     Size: 1-by-1 (scalar)
%                     Units: m/s (meters per second)
% 
%             wstar - Convective velocity scale.
%                     Size: 1-by-1 (scalar)
%                     Units: m/s (meters per second)
% 
%     BLmixingDepth - Boundary layer mixing depth, a key parameter.
%                     Size: 1-by-1 (scalar)
%                     Units: m (meters)
% 
%                 L - Obukhov length.
%                     Size: 1-by-1 (scalar)
%                     Units: m (meters)
% 
%    Outputs:
% 
%           epsilon - The rate of dissipation of turbulent kinetic energy
%                     within the convective boundary layer (CBL).
%                     Size: 1-by-1 (scalar)
%                     Units: m2/s3 (square meters per cubic second)
% 

%% Checks
% No checks

%% Computation
% Calculate the rate of dissipation of TKE
epsilon = 0.4*wstar^3/BLmixingDepth + ustar^3*(1 - z/BLmixingDepth)*(1 - 15*z/L)^-0.25/(0.4*z);