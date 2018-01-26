function [x, y, z] = enu2ecef (e, n, u, lat0, lon0, h0, ell, ...
                               angleut)
                          
  if nargin<7, ell = get_ellipsoid(); end
  if nargin<8, angleut = 'd'; end

  [x0, y0, z0] = geodetic2ecef(lat0, lon0, h0, ell, angleut);
  [dx, dy, dz] = enu2uvw(e, n, u, lat0, lon0, angleut);
  
   x = x0 + dx;
   y = y0 + dy;
   z = z0 + dz;
end