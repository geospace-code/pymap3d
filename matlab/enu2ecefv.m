function [u, v, w] = enu2ecefv(uEast, vNorth, wUp, lat0, lon0, angleut)

  if nargin < 6 
      angleut = 'd';
  end
      
  if strcmpi(angleut(1),'r')
    lat0 = rad2deg(lat0);
    lon0 = rad2deg(lon0);
  end

  t = cos(lat0) .* wUp - sin(lat0) .* vNorth;
  w = sin(lat0) .* wUp + cos(lat0) .* vNorth;

  u = cos(lon0) .* t - sin(lon0) .* uEast;
  v = sin(lon0) .* t + cos(lon0) .* uEast;
  
end