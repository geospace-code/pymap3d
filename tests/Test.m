% Copyright (c) 2014-2018 Michael Hirsch, Ph.D.
%
% Redistribution and use in source and binary forms, with or without modification, are permitted provided that the following conditions are met:
% 1. Redistributions of source code must retain the above copyright notice, this list of conditions and the following disclaimer.
% 2. Redistributions in binary form must reproduce the above copyright notice, this list of conditions and the following disclaimer in the documentation and/or other materials provided with the distribution.
% THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

function Test()

fpath = fileparts(mfilename('fullpath'));
addpath([fpath,filesep,'../matlab'])

E = wgs84Ellipsoid();
%% reference inputs
az = 33; el=70; srange = 1e3;
lat = 42; lon= -82; alt = 200;
%% reference outputs
er = 186.277521; nr = 286.84222; ur = 939.69262; % aer2enu
xl = 660.930e3; yl = -4701.424e3; zl = 4246.579e3; % aer2ecef
lat1 = 42.0026; lon1 = -81.9978; alt1 = 1.1397e3; % aer2geodetic
x0 = 660.675e3; y0 = -4700.949e3; z0 = 4245.738e3; % geodetic2ecef, ecef2geodetic
%% tests

%% aer2ecef contains:
[x1,y1,z1] = geodetic2ecef(E,lat,lon,alt);
assert_allclose([x1,y1,z1],[x0,y0,z0])

[e1,n1,u1] = aer2enu(az, el, srange);
assert_allclose([e1,n1,u1], [er,nr,ur])

[x2,y2,z2] = aer2ecef(az,el,srange,lat,lon,alt,E);
assert_allclose([x2,y2,z2], [xl,yl,zl])

%% ecef2geodetic is self-contained, iterative algorithm.
[lat2, lon2, alt2] = ecef2geodetic(E, x1, y1, z1); % round-trip
assert_allclose([lat2, lon2, alt2], [lat, lon, alt])

[az2, el2, rng2] = enu2aer(e1,n1,u1); % round-trip
assert_allclose([az2,el2,rng2],[az,el,srange])

[az3, el3, rng3] = ecef2aer(x2,y2,z2, lat,lon,alt, E); % round-trip
assert_allclose([az3,el3,rng3], [az,el,srange])


[lat3,lon3,alt3] = aer2geodetic(az,el,srange,lat,lon,alt, E);
assert_allclose([lat3,lon3,alt3], [lat1, lon1, alt1])

[e2, n2, u2] = geodetic2enu(lat3, lon3, alt3, lat, lon, alt, E);
assert_allclose([e2,n2,u2],[e1,n1,u1])

[az4, el4, rng4] = geodetic2aer(lat3,lon3,alt3,lat,lon,alt, E); % round-trip
assert_allclose([az4,el4,rng4], [az,el,srange])
%%
[x3, y3, z3] = enu2ecef(e1,n1,u1,lat,lon,alt, E);
assert_allclose([x3,y3,z3],[x2,y2,z2])

[lat4, lon4, alt4] = enu2geodetic(e2,n2,u2,lat,lon,alt, E); % round-trip
assert_allclose([lat4, lon4, alt4],[lat3, lon3, alt3])

[e3,n3,u3] = ecef2enu(x3,y3,z3,lat,lon,alt, E); % round-trip
assert_allclose([e3,n3,u3],[e1,n1,u1])

end % function


function rad = deg2rad(deg) %#ok<DEFNU>
% for Octave < 4.0
  rad = deg * (pi / 180);

end
