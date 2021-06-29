function [t, x, y, z, up, vp, wp] = LSstep(t, x, y, z, up, vp, wp, station)
% 
% Matt Werner (m.werner@vt.edu) - June 24, 2021
% 
% Evaluate a single time step for a particle flowing through the atmosphere
% experiencing wind velocity components u, v, and w in the respective axis
% directions x, y, and z. The wind velocities (u, v, w) and temperature T
% are to be interpreted as being evaluated at the position (x, y, z) at
% time t.
% 
%    Inputs:
% 
%                 t - Current time relative to when the simulation began.
%                     Size: 1-by-1 (scalar)
%                     Units: s (seconds)
% 
%           x, y, z - Current position of the fluid particle at time t.
%                     Size: 1-by-1 (scalar)
%                     Units: m (meters)
% 
%        up, vp, wp - Current velocity fluctuations of the fluid particle
%                     at time t.
%                     Size: 1-by-1 (scalar)
%                     Units: m/s (meters per second)
% 
%           station - Post-processed data from the measurement system
%                     regarding statistics for the wind (as a function of
%                     height z) including:
%                     1. Elevations (z) for which measurements were taken
%                     2. First elevation at which a NaN value appears
%                     3. Horizontal (u, v) and vertical (w) wind velocities
%                     4. The respective mean values (U, V, W) and
%                        fluctuations (u', v', w')
%                     5. The respective standard deviations (in u, v, w)
%                     6. The respective variances (squared standard dev.)
%                     7. Time-average velocity covariances (u'v' and u'w')
%                     Note: The quantities are to be evaluated with `ppval'
%                           for ease in implementing and choosing the
%                           interpolation method as well as evaluating
%                           derivatives.
%                     Size: 1-by-1 (structure)
%                     Units: ? (SI)
%                     Extended description of fields:
%                           z - Elevations at which measurements are taken
%                               at the measurement station. Most of the
%                               remaining fields will be a function of this
%                               vector (in the form of interp1(..., 'pp')).
%                               Note: Some quantities will fail to be
%                                     measured at generally different
%                                     elevations and above. Thus, it's
%                                     beneficial to obtain a clean set of
%                                     data for most unrestricted
%                                     propagation.
%                               Size: n-by-1 (vector)
%                               Units: m (meters)
%                        zNaN - The lowest elevation at which any measured
%                               value below failed. A simulation of a
%                               particle exceeding this elevation in
%                               simulation should stop.
%                               Size: 1-by-1 (scalar)
%                               Units: m (meters)
%                     u, v, w - Total wind velocity at the measurement
%                               station. For a wind field assumed to be
%                               homogeneous in the downrange direction and
%                               time (that is, u = u(z), v = v(z), and 
%                               w = w(z), i.e. independent of x, y, and t),
%                               these quantities correspond to the wind
%                               field everywhere at the measurement time
%                               t0.
%                               Size: 1-by-1 (structure)
%                               Units: m/s (meters per second)
%                     U, V, W - Mean quantities corresponding to (u, v, w).
%                               The average is performed over the
%                               elevations z at which each were measured.
%                               As such, these quantities are constant
%                               scalars.
%                               Size: 1-by-1 (scalar)
%                               Units: m/s (meters per second)
%                  up, vp, wp - Fluctuating quantities corresponding to 
%                               (u, v, w). These quantities are simply the
%                               difference between the wind velocities and
%                               mean values, though are provided anyway for
%                               clarity since memory is not a concern.
%                               Size: 1-by-1 (structure)
%                               Units: m/s (meters per second)
%            STDu, STDv, STDw - Standard deviations in the measurements of
%                               (u, v, w) at the measurement station taken
%                               at the measurement time t0.
%                               Size: 1-by-1 (structure)
%                               Units: m/s (meters per second)
%            VARu, VARv, VARw - Variance in the measurements of (u, v, w)
%                               at the measurement station taken at the
%                               measurement time t0. These quantities are
%                               simply the standard deviation squared, but
%                               (like the fluctuations) are provided anyway
%                               for clarity since memory is not a concern.
%                               Size: 1-by-1 (structure)
%                               Units: m2/s2 (square meter per square second)
%                  upvp, upwp - Time-averaged Eulerian velocity
%                               covariances.
%                               Size: 1-by-1 (structure)
%                               Units: m2/s2 (square meter per square second)
%   dVARudz, dVARvdz, dVARwdz - Derivative of the variance (in z) in the
%                               measurements of (u, v, w) at the
%                               measurement station taken at the
%                               measurement time t0.
%                               Size: 1-by-1 (structure)
%                               Units: m/s2 (square meter per square second)
%            dupvpdz, dupwpdz - Derivative of time-averaged Eulerian velocity
%                               covariances (in z).
%                               Size: 1-by-1 (structure)
%                               Units: m2/s2 (square meter per square second)
% 
%    Outputs:
% 
%                 t - Future time relative to when the simulation began.
%                     This quantity is equivalent to t + dt, where t here
%                     is the input t.
%                     Size: 1-by-1 (scalar)
%                     Units: s (seconds)
% 
%           x, y, z - Future position of the fluid particle at time t + dt.
%                     Size: 1-by-1 (scalar)
%                     Units: m (meters)
% 
%        up, vp, wp - Future velocity fluctuations of the fluid particle
%                     at time t + dt.
%                     Size: 1-by-1 (scalar)
%                     Units: m/s (meters per second)
% 

%% No checks

%% Extract necessary values at the current elevation z
% Get the mean wind velocity (constants)
U = station.U;
V = station.V;
W = station.W;

% Interpolate the standard deviation
STDu = interp1(station.z, station.STDu, z, 'spline', 'extrap');
STDv = interp1(station.z, station.STDv, z, 'spline', 'extrap');
STDw = interp1(station.z, station.STDw, z, 'spline', 'extrap');

% Interpolate the variance
VARu = STDu^2;
VARv = STDv^2;
VARw = STDw^2;

% Interpolate the covariance
meanupwp = interp1(station.z, station.upwp, z, 'spline', 'extrap');

% Interpolate the derivative of variance
dVARudz = ppval(station.dVARudz, z);
dVARvdz = ppval(station.dVARvdz, z);
dVARwdz = ppval(station.dVARwdz, z);

% Interpolate the derivative of covariance
dmeanupwpdz = ppval(station.dupwpdz, z);

%% Computation
tau = 1; % TEMPORARY % TEMPORARY % TEMPORARY % TEMPORARY % TEMPORARY
% Compute the Langevin coefficients used for propagation for each of the
% three directions (x, y, z)
A = 2*(VARu*VARw - meanupwp^2);
% x direction
bu = 2*VARw/tau;
au = (bu^2/A)*(meanupwp*wp - VARw*up) + 0.5*dmeanupwpdz + ...
        (1/A)*(VARw*dVARudz*up*wp - meanupwp*dVARudz*wp^2 - ...
               meanupwp*dmeanupwpdz*up*wp + VARu*dmeanupwpdz*wp^2);
           
% y direction
bv = bu;
av = -0.5*bv^2*vp/VARv + 0.5*dVARvdz*vp*wp/VARv;

% z direction
bw = bu;
aw = (bw^2/A)*(meanupwp*up - VARu*wp) + 0.5*dVARwdz + ...
        (1/A)*(VARw*dmeanupwpdz*up*wp - meanupwp*dmeanupwpdz*wp^2 - ...
               meanupwp*dVARwdz*up*wp + VARu*dVARwdz*wp^2);
           
%% Dynamics
% Calculate the time step and estimate the change in x iteratively until
% the Courant condition,
%                       dt < dx / (U + 3*STDu),
% is satisfied. This time step will then be used to calculate the remainder
% of the dynamics.
dt = 0.01*tau;
while true
    % Create a random 'infinitesimal' with mean 0 and variance dt. This
    % quantity is created with a new random number upon each iteration
    dxiu = randn*sqrt(dt);
    % Write the stochastic dynamics for the velocity fluctuation in the x
    % direction
    dup = au*dt + bu*dxiu;
    % Write the dynamics for the particle's position in the x direction
    % WITHOUT assigning up = up + dup since the loop will generally take
    % more than 1 iteration to break (don't want accumulation of failed
    % attempts)
    dx = (up + dup + U)*dt;
    
    % Check the Courant condition
    if (dt < dx / (U + 3*STDu))
        % Leave the loop with this time step
        break
    end
    
    % Decrement the time step to obtain finer resolution in the flow
    % characteristics
    dt = 0.5*dt;
    
    % Check if the time step is finite (nonzero)
    if (dt < eps)
        error("Time step decreased to 0.")
    end
end

% Continue evaluating the other dynamics with this time step assuming that
% the coordinate frame is aligned with the mean direction of the wind (i.e.
% V = 0)
% 
% y direction
dxiv = randn*sqrt(dt);
dvp = av*dt + bv*dxiv;
dy = vp*dt;
% z direction
dxiw = randn*sqrt(dt);
dwp = aw*dt + bw*dxiw;
dz = wp*dt;

% Increment the values t, (x, y, z), and (u', v', w') due to this time step
t = t + dt;
x = x + dx;
y = y + dy;
z = z + dz;
up = up + dup;
vp = vp + dvp;
wp = wp + dwp;
