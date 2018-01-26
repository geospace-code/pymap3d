clear
fpath = fileparts(mfilename('fullpath'));
addpath([fpath,filesep,'../matlab'])
ell = get_ellipsoid();

%% reference inputs
az = 33; el=70; srange = 1000;
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

[e,n,u] = aer2enu(az, el, srange);
assert_allclose([e,n,u], [a2e,a2n,a2u])

[ev,nv,uv] = ecef2enuv(vx,vy,vz,lat,lon);
assert_allclose([ev,nv,uv],[ve,vn,vu])

[x,y,z] = aer2ecef(az,el,srange,lat,lon,alt);
assert_allclose([x,y,z], [a2x,a2y,a2z])

%% ecef2geodetic is self-contained, iterative algorithm.
[lat2, lon2, alt2] = ecef2geodetic(x1, y1, z1);
%[lat2, lon2, alt2] = ecef2geodetic(g2x, g2y, g2z);
assert_allclose([lat2, lon2, alt2], [lat, lon, alt])

[lat3,lon3,alt3] = aer2geodetic(az,el,srange,lat,lon,alt);
assert_allclose([lat3,lon3,alt3], [a2la, a2lo, a2a], [], 0.01)
return
[azt, elt, rngt] = ecef2aer(x,y,z,lat,lon,alt,ell,'degrees');
fprintf('\necef2aer %f %f %f',azt,elt,rngt)




[g2az, g2el,g2r] = geodetic2aer(lat2,lon2,alt2,lat,lon,alt,ell,'degrees');
fprintf('geodetic2aer %f %f %f\n',g2az,g2el,g2r)

[ee,en,eu] = ecef2enu(x,y,z,lat,lon,alt,ell);
fprintf('ecef2enu %f %f %f\n',ee,en,eu)

[e2az, e2el, e2rn] = enu2aer(e,n,u,'degrees');
fprintf('enu2aer %f %f %f\n',e2az,e2el,e2rn)

[e2la, e2lo, e2al] = enu2geodetic(e,n,u,lat,lon,alt,ell);
fprintf('enu2geodetic %f %f %f\n',e2la,e2lo,e2al)

[e2x, e2y, e2z ] = enu2ecef(x,y,z,lat,lon,alt,ell,'degrees')

fprintf('\n')
