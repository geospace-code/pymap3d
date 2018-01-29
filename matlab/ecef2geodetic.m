function [lat,lon,alt] = ecef2geodetic(spheroid, x, y, z,  angleUnit)
% function [lat,lon,alt] = ecef2geodetic(spheroid, x, y, z, angleUnit)
%
% Inputs
% ------
% x,y,z:  ECEF coordinates of test point(s) (meters)
% spheroid: referenceEllipsoid parameter struct
% angleUnit: string for angular units. Default 'd': degrees
%
% Outputs
% -------
% lat,lon, alt:  ellipsoid geodetic coordinates of point(s) (degrees, degrees, meters)
%
% also see: http://www.oc.nps.edu/oc2902w/coord/coordcvt.pdf
% Fortran reference at bottom of: http://www.astro.uni.torun.pl/~kb/Papers/geod/Geod-BG.htm

  if isempty(spheroid), spheroid = wgs84Ellipsoid(); end

  % Algorithm is based on 
  % http://www.astro.uni.torun.pl/~kb/Papers/geod/Geod-BG.htm
  % This algorithm provides a converging solution to the latitude equation
  % in terms of the parametric or reduced latitude form (v)
  % This algorithm provides a uniform solution over all latitudes as it does
  % not involve division by cos(phi) or sin(phi)
   a = spheroid.SemimajorAxis; 
   b = spheroid.SemiminorAxis;
  r = hypot(x, y);
  % Constant required for Latitude equation
  rho = atan2(b * z, a * r);  
  % Constant required for latitude equation
  c = (a^2 - b^2) ./ hypot(a*r, b*z);  
   count = 0;
% Starter for the Newtons Iteration Method
  vnew = atan2(a * z, b * r);  
% Initializing the parametric latitude
   v = 0;
  while v ~= vnew && count < 5
 %   disp([v,vnew])
     v = vnew;
    % Derivative of latitude equation
    w = 2 * (cos (v - rho) - c .* cos(2 * v)); 
    % Newtons Method for computing iterations  
    vnew = v - ((2 * sin (v - rho) - c .* sin(2 * v)) ./ w);  
    count = count+1;
  end %while

%% Computing latitude from the root of the latitude equation
  lat = atan2(a * tan (vnew), b); 
  % Computing longitude
  lon = atan2(y, x); 
 % Computing h from latitude obtained 
  alt = ((r - a * cos (vnew)) .* cos (lat)) +  ...
      ((z - b * sin (vnew)) .* sin (lat));
      
  if nargin < 5 || isempty(angleUnit) || strcmpi(angleUnit(1),'d')
    lat = rad2deg(lat);
    lon = rad2deg(lon);
  end

 end % function

% Copyright (c) 2014-2018 Michael Hirsch, Ph.D.
% Copyright (c) 2013, Felipe Geremia Nievinski
%
% Redistribution and use in source and binary forms, with or without modification, are permitted provided that the following conditions are met:
% 1. Redistributions of source code must retain the above copyright notice, this list of conditions and the following disclaimer.
% 2. Redistributions in binary form must reproduce the above copyright notice, this list of conditions and the following disclaimer in the documentation and/or other materials provided with the distribution.
% THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
