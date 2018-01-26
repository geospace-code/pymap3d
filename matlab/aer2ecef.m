function [x,y,z] = aer2ecef(az, el, slantRange, lat0, lon0, h0, ell, ...
                             angleut)
                             
  if nargin < 7 || isempty(ell), ell=get_ellipsoid(); end
  if nargin < 8, angleut = 'd'; end

  %% Origin of the local system in geocentric coordinates.
  [x0, y0, z0] = geodetic2ecef(lat0, lon0, h0, ell, angleut);
  %% Convert Local Spherical AER to ENU
  [e, n, u] = aer2enu(az, el, slantRange, angleut);
  %% Rotating ENU to ECEF
  [dx, dy, dz] = enu2ecefv(e, n, u, lat0, lon0, angleut);
  %% Origin + offset from origin equals position in ECEF
  x = x0 + dx;
  y = y0 + dy;
  z = z0 + dz;

end
