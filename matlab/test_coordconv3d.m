% cross-verify with matlab

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
ell = referenceEllipsoid('wgs84');
[et,nt,ut] = aer2enu(az,el,srange);
fprintf('aer2enu %f %f %f ',et,nt,ut)

[xt,yt,zt] = aer2ecef(az,el,srange,lat,lon,alt,ell);
fprintf('\naer2ecef %f %f %f ',xt,yt,zt)

[latt,lont,altt] = aer2geodetic(az,el,srange,lat,lon,alt,ell,'degrees');
fprintf('\naer2geodetic %f %f %f ',latt,lont,altt)

[azt, elt, rngt] = ecef2aer(x,y,z,lat,lon,alt,ell,'degrees');
fprintf('\necef2aer %f %f %f',azt,elt,rngt)


[g2x,g2y,g2z] = geodetic2ecef(deg2rad(lat),deg2rad(lon),alt,ell);
fprintf('\ngeodetic2ecef %f %f %f',g2x,g2y,g2z)

[g2az, g2el,g2r] = geodetic2aer(lat2,lon2,alt2,lat,lon,alt,ell,'degrees');
fprintf('\ngeodetic2aer %f %f %f',g2az,g2el,g2r)

[ee,en,eu] = ecef2enu(x,y,z,lat,lon,alt,ell);
fprintf('\necef2enu %f %f %f',ee,en,eu)

[ec2az, ec2el, ec2rn] = ecef2geodetic(x,y,z,ell);
fprintf('\necef2geodetic %f %f %f',rad2deg(ec2az),rad2deg(ec2el),ec2rn)

[e2az, e2el, e2rn] = enu2aer(e,n,u,'degrees');
fprintf('\nenu2aer %f %f %f',e2az,e2el,e2rn)

[e2la, e2lo, e2al] = enu2geodetic(e,n,u,lat,lon,alt,ell);
fprintf('\nenu2geodetic %f %f %f',e2la,e2lo,e2al)

[e2x, e2y, e2z ] = enu2ecef(x,y,z,lat,lon,alt,ell,'degrees')

fprintf('\n')