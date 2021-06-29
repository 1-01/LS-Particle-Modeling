% 
% Matt Werner (m.werner@vt.edu) - June 28, 2021
% 
% Make a mock station to test the LS model
% 

% Station measurement elevations
station.z = (200:-10:10)';
% First elevation at which a -9999 (NaN) value appears
station.zNaN = NaN;

% Give the horizontal wind speeds (m/s) and directions (deg)
SPEED = [389,365,339,296,280,281,291,356,396,389,362,339,327,312,293,254,219,181,139,125]'/100;
DIR   = [218,219,220,218,222,224,223,230,226,224,224,224,225,226,226,228,228,224,215,213]';

% Place measured wind speed and directions into North and East components
NORTH = SPEED.*cosd(DIR+180);
EAST = SPEED.*sind(DIR+180);

% Determine the mean direction of wind (determining downrange and
% crossrange directions. The downrange direction is the angle MEANDIR
% measured clockwise from North.
MEANDIR = mean(DIR);
% Rotate North and East components into downrange and crossrange
[DOWNWIND, CROSSWIND] = rotateWind(NORTH, EAST, MEANDIR);

% Assign to station
station.u = DOWNWIND;
station.v = CROSSWIND;
station.w = -[7,7,8,8,9,9,8,9,9,7,8,6,0,2,4,7,12,6,3,7]'/100; % /100 for (m/s)

station.U = mean(station.u);
station.V = mean(station.v);
station.W = mean(station.w);

station.up = station.u - station.U;
station.vp = station.v - station.V;
station.wp = station.w - station.W;

% Can't rotate standard deviation like a vector (?) also unsure how to
% compute it, so for now, assume that the standard deviation in N, E, V is
% the same as downwind, crosswind, vertical
station.STDu = [0,0,1,1,1,4,6,4,2,2,1,1,8,6,10,7,10,3,2,0]'/100;
station.STDv = [188,175,162,148,135,122,105,93,100,109,137,132,118,112,104,94,74,65,55,50]'/100;
station.STDw = [20,19,18,17,16,16,18,15,12,13,10,9,9,14,15,13,8,4,2,5]'/100;

station.VARu = station.STDu.^2;
station.VARv = station.STDv.^2;
station.VARw = station.STDw.^2;

% Unsure if U'V' and U'W' in the file are North/East/Vertical or
% downrange/crossrange/vertical and unsure if they are rotatable, so for
% now just assume that they're downrange/crossrange/vertical
station.upvp = [1,0,0,0,4,2,2,2,1,0,2,2,5,5,11,8,11,6,3,1]'/100;
station.upwp = -[135,128,122,115,109,102,81,95,85,132,138,126,114,111,106,100,85,70,39,32]'/100;

% Make a structure to evaluate the derivative of the variances
% u
varustruct = interp1(station.z, station.VARu, 'spline', 'pp');
station.dVARudz = fnder(varustruct,1);
% v
varvstruct = interp1(station.z, station.VARv, 'spline', 'pp');
station.dVARvdz = fnder(varvstruct,1);
% w
varwstruct = interp1(station.z, station.VARw, 'spline', 'pp');
station.dVARwdz = fnder(varwstruct,1);
% u'w'
varupwpstruct = interp1(station.z, station.upwp, 'spline', 'pp');
station.dupwpdz = fnder(varupwpstruct, 1);