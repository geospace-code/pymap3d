function [uEast, vNorth, wUp] = ecef2enuv (u, v, w, lat0, lon0, angleut)

  if nargin<6 || strcmpi(angleut(1), 'd')
    lat0 = deg2rad(lat0);
    lon0 = deg2rad(lon0);
  end

  t      =  cos(lon0) .* u + sin(lon0) .* v;
  uEast  = -sin(lon0) .* u + cos(lon0) .* v;
  wUp    =  cos(lat0) .* t + sin(lat0) .* w;
  vNorth = -sin(lat0) .* t + cos(lat0) .* w;
end
