function [x,y,z] = geodetic2ecef(spheroid, lat, lon, alt,  angleUnit)
% [x,y,z] = geodetic2ecef(spheroid, lat, lon, alt,  angleUnit)
%
% Inputs
% ------
% lat,lon, alt:  ellipsoid geodetic coordinates of point(s) (degrees, degrees, meters)
% spheroid: referenceEllipsoid parameter struct
% angleUnit: string for angular units. Default 'd': degrees
%
% outputs
% -------
% x,y,z:  ECEF coordinates of test point(s) (meters)
%

  if isempty(spheroid), spheroid = wgs84Ellipsoid(); end
  
  if nargin < 5 || isempty(angleUnit) || strcmpi(angleUnit(1),'d')
    lat = deg2rad(lat);
    lon = deg2rad(lon);
  end
 
  
%% Radius of curvature of the prime vertical section
  N = get_radius_normal(lat, spheroid);
%% Compute cartesian (geocentric) coordinates given  (curvilinear) geodetic coordinates.
  
  x = (N + alt) .* cos(lat) .* cos(lon);
  y = (N + alt) .* cos(lat) .* sin(lon);
  z = (N .* (spheroid.SemiminorAxis / spheroid.SemimajorAxis)^2 + alt) .* sin(lat);
        
end
 
% Copyright (c) 2014-2018 Michael Hirsch, Ph.D.
% Copyright (c) 2013, Felipe Geremia Nievinski
%
% Redistribution and use in source and binary forms, with or without modification, are permitted provided that the following conditions are met:
% 1. Redistributions of source code must retain the above copyright notice, this list of conditions and the following disclaimer.
% 2. Redistributions in binary form must reproduce the above copyright notice, this list of conditions and the following disclaimer in the documentation and/or other materials provided with the distribution.
% THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
