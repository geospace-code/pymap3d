function [lat, lon, h] = enu2geodetic (e, n, u, lat0, lon0, ...
                                       h0, ell, angleut)
             
  if nargin<7, ell = get_ellipsoid(); end
  if nargin<8, angleut='d'; end

  [x, y, z] = enu2ecef(e, n, u, lat0, lon0, h0, ell, angleut);
  [lat, lon, h] = ecef2geodetic(x, y, z, ell, angleut);

end