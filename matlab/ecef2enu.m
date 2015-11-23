%%
%This is adapted from Octave Mapping Toolbox by Michael Hirsch, so as to be Matlab and Octave compatible
% Copyright (C) 2013 Felipe G. Nievinski
% Copyright (C) 2013 Sandeep V. Manthi

function [xEast,yNorth,zUp] = ecef2enu (x, y, z, lat0, lon0, h0, ell, angleut)
  if nargin < 7 || isempty (ell)
    ell = get_ellipsoid('wgs84');
  elseif ~isstruct (ell)
    ell = get_ellipsoid (ell);
  endif

  if nargin<8
    angleut='degree';
  end

  [x0, y0, z0] = geodetic2ecef (lat0, lon0, h0, ell, angleut);
  [xEast, yNorth, zUp] = ecef2enuv (x - x0, y - y0, z - z0, lat0, lon0, angleut);
end
