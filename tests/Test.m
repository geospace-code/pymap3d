matlabmap = true;

if ~matlabmap
  fpath = fileparts(mfilename('fullpath'));
  addpath([fpath,filesep,'../matlab'])
  ell = get_ellipsoid();
else
  ell = referenceEllipsoid('wgs84');
end
%% reference inputs
az = 33; el=70; srange = 1000; % aer2enu
lat = 42; lon= -82; alt= 200;
%% reference outputs
a2e = 186.277521; a2n = 286.842228; a2u = 939.692621; % aer2enu
a2x = 660930.192761; a2y= -4701424.222957; a2z= 4246579.604633; % aer2ecef
ec2la= 42.100000; ec2lo= -81.900000; ec2a = 300.000000; %ecef2geodetic
%% tests
[e,n,u] = aer2enu(az, el, srange);
assert_allclose([e,n,u], [a2e,a2n,a2u], [], 0.01)

[x,y,z] = aer2ecef(az,el,srange,lat,lon,alt,ell);
assert_allclose([x,y,z], [a2x,a2y,a2z], [], 0.01)

[lat2, lon2, alt2] = ecef2geodetic(x, y, z,ell);
assert_allclose([lat2, lon2, alt2], [lat, lon, alt])

[lat3,lon3,alt3] = aer2geodetic(az,el,srange,lat,lon,alt,ell);
assert_allclose([lat3,lon3,alt3], [lat, lon, alt], [], 0.01)
return
[azt, elt, rngt] = ecef2aer(x,y,z,lat,lon,alt,ell,'degrees');
fprintf('\necef2aer %f %f %f',azt,elt,rngt)


[g2x,g2y,g2z] = geodetic2ecef(deg2rad(lat),deg2rad(lon),alt,ell);
fprintf('\ngeodetic2ecef %f %f %f',g2x,g2y,g2z)

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
