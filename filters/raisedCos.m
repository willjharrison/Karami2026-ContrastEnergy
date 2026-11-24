function filter = raisedCos(radDist,ctr,cut,logOrLinear)
% filter = raisedCos(radDist,ctr,cut,logOrLinear)
%
% Creates a raised log (or linear) cosine filter, similar to that described
% by Peli (1990), and used by Chung et al. radDist is a matrix of radial
% distances relative to the centre of an image (or 2d fourier transform).
% ctr is the peak frequency of the filter.

if nargin < 4
    logOrLinear = 'log';
end

if strcmp(logOrLinear,'log')
    filter = (1 + cos(pi*((log(radDist)-log(ctr))/(log(cut)-log(ctr)))))/2;
    filter(radDist < cut | radDist > ctr^2/cut) = 0;
else
    
    % Linear-scale raised-cosine filter
    filter = (1 + cos(pi*(rad2deg(wrap(deg2rad(radDist*2-ctr*2))/(cut*2-ctr*2)))))/2;
    filter(abs(rad2deg(wrap(deg2rad(radDist-ctr)))) > abs(cut - ctr)) = 0;

%     filter = (1 + cos(pi*(rad2deg((deg2rad(radDist*2-ctr*2))/(cut*2-ctr*2)))))/2;
%     filter(abs(rad2deg((deg2rad(radDist-ctr)))) > abs(cut - ctr)) = 0;

end


function X = wrap(Y, bound)
% X = WRAP(Y)
%   Maps Y values onto circular space (-PI <= X < PI).
%   X = WRAP(Y, BOUND) specifies alternative bounds (default is PI).
%
%   --> www.paulbays.com

if nargin<2, bound = pi; end

X = mod(Y + bound, bound*2) - bound;
