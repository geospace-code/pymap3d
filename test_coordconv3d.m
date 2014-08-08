% cross-verify with matlab

lat = 42; lon = -82; alt = 200;
az = 45; el = 80; srange = 1000;
 x = 23893;  y=2193; z =  2389239;

%% 
ell = referenceEllipsoid('wgs84');
[et,nt,ut] = aer2enu(az,el,srange);
fprintf('aer2enu %f %f %f ',et,nt,ut)

[xt,yt,zt] = aer2ecef(az,el,srange,lat,lon,alt,ell);
fprintf('\naer2ecef %f %f %f ',xt,yt,zt)

[latt,lont,altt] = aer2geodetic(az,el,srange,lat,lon,alt,ell);
fprintf('\naer2geodetic %f %f %f ',latt,lont,altt)

[azt, elt, rngt] = ecef2aer(x,y,z,lat,lon,alt,ell);
fprintf('\necef2aer %f %f %f',azt,elt,rngt)
display(' ')