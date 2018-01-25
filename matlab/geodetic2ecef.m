%%
%This is adapted from Octave Mapping Toolbox by Michael Hirsch, so as to be Matlab and Octave compatible
% Copyright (C) 2013 Felipe G. Nievinski
% Copyright (C) 2013 Sandeep V. Manthi

function [x,y,z] = geodetic2ecef(phi, lambda, h, ell, angleut)

  % CONVERT_TO_CARTESIAN: Convert to global Cartesian coordinates (ECEF),
  % given geodetic curvilinear coordinates.
  coord_geod = [phi,lambda, h];
  if nargin < 4 || isempty(ell) 
    ell = get_ellipsoid();
  elseif ~isstruct(ell)
    ell = get_ellipsoid(ell);
  end

  if nargin<5
    angleut='degree';
  end

  % Radius of curvature of the prime vertical section
  N = get_radius_normal(coord_geod(:, 1), ell);

  % Some shortnames for variables used often.
  a = ell.a;  
  b = ell.b;
  lat = coord_geod(:, 1);
  lon = coord_geod(:, 2);
  h   = coord_geod(:, 3);

  % Compute cartesian (geocentric) coordinates given 
  % (curvilinear) geodetic coordinates.
  if ~strcmpi(angleut(1),'d')
    lat = rad2deg(lat);
    lon = rad2deg(lon);
  end
  
  x = (N + h) .* cosd(lat) .* cosd(lon);
  y = (N + h) .* cosd(lat) .* sind(lon);
  z = (N .* (b / a)^2 + h) .* sind(lat);
        
 end
