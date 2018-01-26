%%
%This is adapted from Octave Mapping Toolbox by Michael Hirsch, so as to be Matlab and Octave compatible
% Copyright (C) 2013 Felipe G. Nievinski
% Copyright (C) 2013 Sandeep V. Manthi

function [x,y,z] = geodetic2ecef(lat, lon, alt, ell, angleut)

  if nargin < 4 || isempty(ell), ell = get_ellipsoid(); end
  if nargin < 5, angleut='d'; end
 
  if strcmpi(angleut(1),'d')
    lat = deg2rad(lat);
    lon = deg2rad(lon);
  end
 
  
%% Radius of curvature of the prime vertical section
  N = get_radius_normal(lat, ell);
%% Compute cartesian (geocentric) coordinates given  (curvilinear) geodetic coordinates.
  
  x = (N + alt) .* cos(lat) .* cos(lon);
  y = (N + alt) .* cos(lat) .* sin(lon);
  z = (N .* (ell.b / ell.a)^2 + alt) .* sin(lat);
        
 end
