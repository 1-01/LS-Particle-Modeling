function TLW = LagrangianTimeScale(ustar, L, cropHeight, z, SIGW)
% 
% Matt Werner (m.werner@vt.edu) - July 14, 2021
% 
% Calculate the Lagrangian timescale of the flow field according to code
% adapted from Aylor in BASIC.
% 
% The Lagrangian timescale is used in determining the local turbulent
% decorrelation time scale for velocity.
% 
%    Inputs:
% 
%             ustar - Friction velocity.
%                     Size: 1-by-1 (scalar)
%                     Units: m/s (meters per second)
% 
%                 L - Obukhov length.
%                     Size: 1-by-1 (scalar)
%                     Units: m (meters)
% 
%                 z - Current elevation of particle.
%                     Size: 1-by-1 (scalar)
%                     Units: m (meters)
% 
%              SIGW - Standard deviation at height z (used only if
%                     z/cropHeight > 1).
%                     Size: 1-by-1 (scalar)
%                     Units: m/s (meters per second)
% 
%        cropHeight - Total height of the crops (from which the particles
%                     are released. This is NOT the height at which the
%                     particles are released, but rather how tall the crops
%                     are).
%                     Size: 1-by-1 (scalar)
%                     Units: m (meters)
% 
%    Outputs:
% 
%                TL - The Lagrangian timescale.
%                     Size: 1-by-1 (scalar)
%                     Units: ? (SI)
% 

%% Checks
% No checks

%% Computations
% Define the factual Lagrangian timescale
TLW_FACT = 0.5;

% Calculate the Lagrangian timescale at 0 altitude according to the sign of
% the Obukhov length.
% 
% Note: In Aylor's notation, L is MoninL (for the Monin-Obukhov length)
if (L > 0)
    % Stable and neutral conditions
    
    % Define the standard deviation of vertical wind at the ground to be as
    % such above the roughness sublayer (RSL)
    SIGW0 = 1.25*ustar;
    
    % Formula of Wilson (1981)
    % Note: If use 0.4 (for TLW_FACT?), then agrees with Raupach
    TLW0 = TLW_FACT*cropHeight/(1 + 5*cropHeight/L)/SIGW0;
    % Correct the high spores to keep them under control
    TLW0 = TLW0/(1 + 15*cropHeight/1000);
else
    % Unstable conditions
    
    % Set SigW above canopy (from Panofsky 1977)
    Disp_Hgt = evalin('base','z(1)');
    SIGW0 = ustar*(2.2 - 6.6*(cropHeight - Disp_Hgt)/L)^0.333;
    
    % Set the Lagrangian timescale in the canopy using formula of Wilson
    % (1981)(conflicts with advice of Raupauch).
    TLW0 = TLW_FACT*cropHeight*(1 - 6*cropHeight/L)^0.25/SIGW0;
    % Adapt to keep high particles under control
    TLW0 = TLW0/(1 + 15*cropHeight/1000);
end

% Calculate the Lagrangian timescale
z_on_h = z/cropHeight;
if (z_on_h < 0.25)
    TLW = TLW0*(0.1 + 3.6*z_on_h);
elseif (z_on_h < 1)
    TLW = TLW0;
else
    if (L > 0)
        TLW = TLW_FACT*z/(1 + 5*z/L)/SIGW;
        TLW = TLW/(1 + 15*z/1000);
    else
        TLW = TLW_FACT*z*(1 - 6*z/L)^0.25/SIGW;
        TLW = TLW/(1 + 15*z/1000);
    end
    % Ensure the Lagrangian timescale isn't less than that at the surface
    if (TLW < TLW0), TLW = TLW0; end
end