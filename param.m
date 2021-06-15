
function [h, d, h0, z0, zw, L, ustar, v_s, k, np, zmin, zmax, tmax, Q_LS, t0, deltstart, sigu_zw, sigw_zw, sigu_h, sigw_h, ex_h, psi_h, Ubar_h, T_Lh, T_L0, beta, x_sense] = param(exp);
%ENVIRONMENTAL PARAMETERS
h = [0.75 0.85 0.95 0.97 0.91 0.95]; % (m) canopy height
h =  h (exp)

h0 = [.53 .47 .53 .52 .44 .42; ... % release height of particles row1 = lower,
      .89 .79 .80 .78 .74 .71]*h;  % row2 = upper
ustar = [.49 .3 .43 .35 .38 .39]/2;   %friction velocity
d = [.7 .7 .7 .7 .7 0.62]*h;   % according to caption, it gives h/d = 0.7 for (a)-(e), 0.62 for (f)
                                % BUT appendix says d/h
d = d(exp)

L = [-100 -40 -100 -40 -10000 -40];    %monin obhukov constant -40 m for exp 6 in Aylor 2001
%z0 = [.11 .08 .08 .065 .11 .1]*d;  %roughness length, taken to be 0.1*d in sixth experiment (also value is described in pg. 48 of Urban transport book)
z0 = [.11 .08 .08 .065 .11 .1]*h; % According to appendix, it's z0/h, not z0/d
 
x_sense = [1 2 2 2 2 2]; % according to caption, at x=2m for (b)-(f)


h0 = h0(:,exp)
ustar = ustar(exp)
L = L(exp)
z0 = z0(exp)
x_sense = x_sense(exp)


%zi = 1000;  %height of mixed layer (assumed to be 1000m)
zw = z0+h;  %height of roughness sublayer, first paragraph of the appendix


%Dp = 30;    %(microns) average diameter of particle (if dp is 400, many fall, some go up)
v_s = .019;  % .019 from 2001 paper for diameter of lycopodium spores // (3*10^-3*Dp^2)/100 (m/s) approximates settling velocity from Aylor's text eq. 3.13
k = 0.4;    %von karman constant

%% MODEL RUN PARAMETERS
np = 100;   %number of particles released at h0
zmin = 0;   %domain lowest height
zmax = 15; %domain max height
tmax = 500; %time to run simulation
Q_LS = 2;   %spores/second
t0 = 0;     %time at which particles start being released
deltstart = 1/Q_LS; %time gap in seconds between each spore released
sigu_zw = ustar*(4+0.6*(-1000/L)^(2/3))^.5;  %horizontal turb variance at upper boundary of roughness sublayer (not height dependent, used for all z > zw)
sigw_zw = 1.25*ustar*(1-3*((zw-d)/L))^(1/3);    %horizontal turb variance at upper boundary of roughness sublayer
sigu_h = 0.85*sigu_zw;
sigw_h = 0.85*sigw_zw;
% sigu_h = 2*ustar;       %from appendix part c Aylor 2001
% sigw_h = 1.25*ustar*(1-3*(h/L))^(1/3);
ex_h = (1-15*(h-d)/L)^0.25;
psi_h = -2*log((1+ex_h)/2) - log((1+ex_h^2)/2) + 2*atan(ex_h) - pi/2;
Ubar_h = ustar/k*(log((h-d)/z0) + psi_h);     %eq A1 evaluated at crop height Aylor 2001
T_Lh = (0.5*h/sigw_h)*(1-6*h/L)^0.25;
T_L0 = T_Lh/(1+15*h/1000);  %from his code "to keep high particles under control"
beta = 2;

end

