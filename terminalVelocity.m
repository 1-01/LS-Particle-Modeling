function vs = terminalVelocity(particleDiameter, particleDensity, ...
                                 airTemperature,      airDensity)
% 
% Matt Werner (m.werner@vt.edu) - July 6, 2021
% 
% Estimate the terminal velocity of a particle (of some diameter and
% density) within the Stokes flow regime. 
% 
% Atmospheric effects propagate into the determination of the terminal
% velocity through the air density and temperature (viscosity).
% 
% Note: The particle is assumed to have a spherical, (or at most nearly
%       spherical) surface.
% 
% Note: When Aylor says "unit density", he means in units of g/cm3
%       (1 g/cm3 = 1000 kg/m3)
% 
%    Inputs:
% 
%  particleDiameter - Particle diameter.
%                     Size: 1-by-1 (scalar)
%                     Units: microns (micrometers)
% 
%   particleDensity - Particle diameter.
%                     Size: 1-by-1 (scalar)
%                     Units: g/cm3 (grams per cubic centimetre)
% 
%    airTemperature - Air temperature.
%                     Size: 1-by-1 (scalar)
%                     Units: C (Celsius)
% 
%        airDensity - Air density.
%                     Size: 1-by-1 (scalar)
%                     Units: kg/m3 (kilogram per cubic meter)
% 
%    Outputs:
% 
%                vs - Terminal velocity of the particle in still air.
%                     Size: 1-by-1 (scalar)
%                     Units: m/s (meters per second)
% 

%% Checks
% Ensure that the particle density is greater than the air density
assert(1e3*particleDensity > airDensity, "Particle density must be greater than air density")

%% Computation
% Convert particle diameter from microns to meters and particle density
% from g/cm3 to kg/m3
particleDiameter = 1e-6*particleDiameter;
particleDensity  = 1e+3*particleDensity;
% Surface gravity (note: independent of latitude)
g = 9.81; % (m/s2)

% Determine atmospheric viscosity
viscosity = dynamicViscosity(273.15 + airTemperature);

% Spherical particles having a diameter smaller than 50 microns follow a
% direct formula to determine the terminal velocity
if (particleDiameter < 50e-6)
    % Terminal velocity (Aylor Eq. 3.10 and 3.11)
    vs = (particleDensity - airDensity)*particleDiameter^2*g/(18*viscosity);
    % Do NOT correct for particles smaller than 1 micron as the correction
    % is nearly insignificant (< about 1%). If performing the correction,
    % then the correction factor (Cc) would be different from unity (1)
    Cc = 1;
    vs = vs*Cc;
    return
end
% Otherwise, create an implicit function for the terminal velocity and
% solve for it numerically (eq. 3.16, 3.17)
% 
% Cross-sectional area of the particle
Ap = (pi/4)*particleDiameter^2;
% Volume of the particle
Vp = (pi/6)*particleDiameter^3;
% Solve for the terminal velocity numerically
Re = @(u) airDensity*u*particleDiameter/viscosity;
cd = @(u) 24./Re(u) .* (1 + 0.158*Re(u).^(2/3));
f = @(vs) vs.^2 - 2*(particleDensity - airDensity)*Vp*g./(airDensity*cd(vs)*Ap);
vs = fzero(f, 0.25);