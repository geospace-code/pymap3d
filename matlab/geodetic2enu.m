function [e, n, u] = geodetic2enu(lat, lon, alt, lat0, lon0, alt0, spheroid, angleUnit)
% function [e, n, u] = geodetic2enu(lat, lon, alt, lat0, lon0, alt0, spheroid, angleUnit)
%
% Inputs
% ------
% lat,lon, alt:  ellipsoid geodetic coordinates of point under test (degrees, degrees, meters)
% lat0, lon0, alt0: ellipsoid geodetic coordinates of observer/reference (degrees, degrees, meters)
% spheroid: referenceEllipsoid parameter struct
% angleUnit: string for angular units. Default 'd': degrees
%
% outputs
% -------
% e,n,u:  East, North, Up coordinates of test points (meters)


  if nargin < 7, spheroid = []; end
  if nargin < 8, angleUnit = []; end
  
  [x1,y1,z1] = geodetic2ecef(spheroid, lat,lon,alt,angleUnit);
  [x2,y2,z2] = geodetic2ecef(spheroid, lat0,lon0,alt0,angleUnit);
  
  dx = x1-x2;
  dy = y1-y2;
  dz = z1-z2;
  
  [e, n, u] = ecef2enuv(dx, dy, dz, lat0, lon0, angleUnit);
  

end % function

% Copyright (c) 2014-2018 Michael Hirsch, Ph.D.
% Copyright (c) 2013, Felipe Geremia Nievinski
%
% Redistribution and use in source and binary forms, with or without modification, are permitted provided that the following conditions are met:
% 1. Redistributions of source code must retain the above copyright notice, this list of conditions and the following disclaimer.
% 2. Redistributions in binary form must reproduce the above copyright notice, this list of conditions and the following disclaimer in the documentation and/or other materials provided with the distribution.
% THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
