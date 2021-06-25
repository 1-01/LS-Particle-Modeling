function meanxy = meanProduct(x, y)
% 
% Matt Werner (m.werner@vt.edu) - June 25, 2021
% 
% Calculate the mean element-wise product of two n-by-1 vectors.
% 
%    Inputs:
% 
%              x, y - Quantities whose mean product is desired to be
%                     calculated. Note that the mean product is calculated
%                     element-wise; as such, the dimensions of x and y must
%                     match.
%                     Size: n-by-1 (vector)
%                     Units: ? (SI)
% 
%    Outputs:
% 
%            meanxy - The mean product of the quantities x and y.
%                     Size: 1-by-1 (scalar)
%                     Units: ? (SI)
% 

%% Checks
% Ensure x and y are of the same dimension
assert(all(size(x) == size(y)), "Dimensions must match.")
assert(size(x, 2) == 1, "Quantities must be column vectors")

%% Computation
% Calculate the mean product
meanxy = x'*y/numel(x);