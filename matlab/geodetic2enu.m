function [xEast, yNorth, zUp] = geodetic2enu(lat, lon, h, lat0, lon0, h0, ell, angleut)

  if nargin<8
    angleut = 'd';
  end
  
  [x1,y1,z1] = geodetic2ecef(lat,lon,h,ell,angleut);
  [x2,y2,z2] = geodetic2ecef(lat0,lon0,h0,ell,angleut);
  dx = x1-x2;
  dy = y1-y2;
  dz = z1-z2;
  [xEast, yNorth, zUp] = ecef2enuv(dx, dy, dz, lat0, lon0, angleut);

endfunction
