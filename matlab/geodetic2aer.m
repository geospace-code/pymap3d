function [az, el, slantRange] = geodetic2aer(lat, lon, h, lat0, lon0, h0, ell, angleut)

  if nargin<8
    angleut = 'd';
  end

  [e, n, u] = geodetic2enu(lat, lon, h, lat0, lon0, h0, ell,angleut );
  [az, el, slantRange] = enu2aer(e, n, u, angleut);
  
end
