clear
fpath = fileparts(mfilename('fullpath'));
addpath([fpath,filesep,'../matlab'])
ell = get_ellipsoid();

%% reference inputs
az = 33; el=70; srange = 1e3;
lat = 42; lon= -82; alt = 200;
vx = 5; vy=3; vz=2;
%% reference outputs
a2e = 186.277521; a2n = 286.84222; a2u = 939.69262; % aer2enu
a2x = 660.930e3; a2y= -4701.424e3; a2z= 4246.579e3; % aer2ecef
a2la = 42.0026; a2lo=-81.9978; a2a=1.1397e3; % aer2geodetic
g2x = 660.675e3; g2y= -4700.949e3; g2z= 4245.738e3; % geodetic2ecef, ecef2geodetic
ve=5.368859646588048; vn= 3.008520763668120; vu= -0.352347711524077; % ecef2enuv
%% tests

%% aer2ecef contains:
[x1,y1,z1] = geodetic2ecef(lat,lon,alt);
assert_allclose([x1,y1,z1],[g2x,g2y,g2z])

[e1,n1,u1] = aer2enu(az, el, srange);
assert_allclose([e1,n1,u1], [a2e,a2n,a2u])

[ev,nv,uv] = ecef2enuv(vx,vy,vz,lat,lon);
assert_allclose([ev,nv,uv],[ve,vn,vu])

[x2,y2,z2] = aer2ecef(az,el,srange,lat,lon,alt);
assert_allclose([x2,y2,z2], [a2x,a2y,a2z])

%% ecef2geodetic is self-contained, iterative algorithm.
[lat2, lon2, alt2] = ecef2geodetic(x1, y1, z1); % round-trip
%[lat2, lon2, alt2] = ecef2geodetic(g2x, g2y, g2z);
assert_allclose([lat2, lon2, alt2], [lat, lon, alt])

[az2, el2, rng2] = enu2aer(e1,n1,u1); % round-trip
assert_allclose([az2,el2,rng2],[az,el,srange])

[az3, el3, rng3] = ecef2aer(x2,y2,z2, lat,lon,alt); % round-trip 
assert_allclose([az3,el3,rng3], [az,el,srange])


[lat3,lon3,alt3] = aer2geodetic(az,el,srange,lat,lon,alt);
assert_allclose([lat3,lon3,alt3], [a2la, a2lo, a2a])

[e2, n2, u2] = geodetic2enu(lat3, lon3, alt3, lat, lon, alt);
assert_allclose([e2,n2,u2],[e1,n1,u1])

[az4, el4, rng4] = geodetic2aer(lat3,lon3,alt3,lat,lon,alt); % round-trip
assert_allclose([az4,el4,rng4], [az,el,srange])
%% 
[x3, y3, z3] = enu2ecef(e1,n1,u1,lat,lon,alt);
assert_allclose([x3,y3,z3],[x2,y2,z2])

[lat4, lon4, alt4] = enu2geodetic(e2,n2,u2,lat,lon,alt); % round-trip
assert_allclose([lat4, lon4, alt4],[lat3, lon3, alt3])

return

[ee,en,eu] = ecef2enu(x,y,z,lat,lon,alt,ell);
fprintf('ecef2enu %f %f %f\n',ee,en,eu)



