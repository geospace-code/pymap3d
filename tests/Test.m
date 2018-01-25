matlabmap = false;

if ~matlabmap
  fpath = fileparts(mfilename('fullpath'));
  addpath([fpath,filesep,'../matlab'])
  ell = get_ellipsoid();
else
  ell = referenceEllipsoid('wgs84');
end

lat = 42; lon = -82; alt = 200;

lat2 = 42.1; lon2 = -81.9; alt2 = 1300;
e =  8.273771039503677e+03;
n =  1.111452002615149e+04;
u =  1.084939260985176e+03;

az = 33; el = 70; srange = 1000;
 x = 6.678411289903646e+05;  
 y =  -4.692496355102768e+06; 
 z = 4.254052899714093e+06;

%%
[et,nt,ut] = aer2enu(az,el,srange);
fprintf('aer2enu %f %f %f \n',et,nt,ut)

[xt,yt,zt] = aer2ecef(az,el,srange,lat,lon,alt,ell);
fprintf('aer2ecef %f %f %f \n',xt,yt,zt)

[latt,lont,altt] = aer2geodetic(az,el,srange,lat,lon,alt,ell,'degrees');
fprintf('aer2geodetic %f %f %f \n',latt,lont,altt)

[azt, elt, rngt] = ecef2aer(x,y,z,lat,lon,alt,ell,'degrees');
fprintf('\necef2aer %f %f %f',azt,elt,rngt)


[g2x,g2y,g2z] = geodetic2ecef(deg2rad(lat),deg2rad(lon),alt,ell);
fprintf('\ngeodetic2ecef %f %f %f',g2x,g2y,g2z)

[g2az, g2el,g2r] = geodetic2aer(lat2,lon2,alt2,lat,lon,alt,ell,'degrees');
fprintf('geodetic2aer %f %f %f\n',g2az,g2el,g2r)

[ee,en,eu] = ecef2enu(x,y,z,lat,lon,alt,ell);
fprintf('ecef2enu %f %f %f\n',ee,en,eu)

[ec2az, ec2el, ec2rn] = ecef2geodetic(x,y,z,ell);
fprintf('ecef2geodetic %f %f %f\n',rad2deg(ec2az),rad2deg(ec2el),ec2rn)

[e2az, e2el, e2rn] = enu2aer(e,n,u,'degrees');
fprintf('enu2aer %f %f %f\n',e2az,e2el,e2rn)

[e2la, e2lo, e2al] = enu2geodetic(e,n,u,lat,lon,alt,ell);
fprintf('enu2geodetic %f %f %f\n',e2la,e2lo,e2al)

[e2x, e2y, e2z ] = enu2ecef(x,y,z,lat,lon,alt,ell,'degrees')

fprintf('\n')
