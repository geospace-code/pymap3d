function [lat, lon, d] = lookAtSpheroid(lat0, lon0, h0, az, tilt, spheroid, angleUnit)
% Calculates line-of-sight intersection with Earth (or other ellipsoid) surface from above surface ./ orbit
%
% Args:
% lat0, lon0: latitude and longitude of starting point
% h0: altitude of starting point in meters
% az: azimuth angle of line-of-sight, clockwise from North
% tilt: tilt angle of line-of-sight with respect to local vertical (nadir = 0)
%
% Returns:
% lat, lon: latitude and longitude where the line-of-sight intersects with the Earth ellipsoid
% d: slant range in meters from the starting point to the intersect point
%
%  Values will be NaN if the line of sight does not intersect.
%
% Algorithm based on https://medium.com/@stephenhartzell/satellite-line-of-sight-intersection-with-earth-d786b4a6a9b6 Stephen Hartzell
narginchk(5,7)
if nargin < 6 || isempty(spheroid), spheroid = wgs84Ellipsoid(); end
if nargin < 7 || isempty(angleUnit), angleUnit='d'; end

validateattributes(lat0, {'numeric'}, {'real','>=',-90,'<=',90})
validateattributes(lon0, {'numeric'}, {'real','finite'})
validateattributes(h0, {'numeric'}, {'real','nonnegative','finite'})
validateattributes(az, {'numeric'}, {'real','finite'})
validateattributes(tilt, {'numeric'}, {'real','nonnegative','<=',180})
validateattributes(spheroid,{'struct'},{'nonempty'})
validateattributes(angleUnit,{'string','char'},{'scalar'})
%% computation

if strcmpi(angleUnit(1),'d')
  el = tilt - 90;
else
  el = tilt - pi/2;
end

a = spheroid.SemimajorAxis;
b = a;
c = spheroid.SemiminorAxis;

[e, n, u] = aer2enu(az, el, 1., angleUnit);  % fixed 1 km slant range
[u, v, w] = enu2uvw(e, n, u, lat0, lon0, angleUnit);
[x, y, z] = geodetic2ecef([], lat0, lon0, h0, angleUnit);

value = -a.^2 .* b.^2 .* w .* z - a.^2 .* c.^2 .* v .* y - b.^2 .* c.^2 .* u .* x;
radical = a.^2 .* b.^2 .* w.^2 + a.^2 .* c.^2 .* v.^2 - a.^2 .* v.^2 .* z.^2 + 2 .* a.^2 .* v .* w .* y .* z - ...
           a.^2 .* w.^2 .* y.^2 + b.^2 .* c.^2 .* u.^2 - b.^2 .* u.^2 .* z.^2 + 2 .* b.^2 .* u .* w .* x .* z - ...
           b.^2 .* w.^2 .* x.^2 - c.^2 .* u.^2 .* y.^2 + 2 .* c.^2 .* u .* v .* x .* y - c.^2 .* v.^2 .* x.^2;

magnitude = a.^2 .* b.^2 .* w.^2 + a.^2 .* c.^2 .* v.^2 + b.^2 .* c.^2 .* u.^2;

% Return nan if radical < 0 or d < 0 because LOS vector does not point towards Earth 
d = (value - a .* b .* c .* sqrt(radical)) ./ magnitude;
d(radical < 0 | d < 0) = nan; % separate line

% altitude should be zero
[lat, lon] = ecef2geodetic([], x + d .* u, y + d .* v, z + d .* w, angleUnit);

end % function


% Copyright (c) 2018 Michael Hirsch, Ph.D.
%
% Redistribution and use in source and binary forms, with or without modification, are permitted provided that the following conditions are met:
% 1. Redistributions of source code must retain the above copyright notice, this list of conditions and the following disclaimer.
% 2. Redistributions in binary form must reproduce the above copyright notice, this list of conditions and the following disclaimer in the documentation and/or other materials provided with the distribution.
% THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE
