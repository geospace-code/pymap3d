function [uEast, vNorth, wUp] = ecef2enuv (u, v, w, lat0, lon0, angleut)
  if nargin<6
    angleut='d';
  end

  if ~strcmpi(angleut(1), 'd')
    lat0 = rad2deg(lat0);
    lon0 = rad2deg(lon0);
  end

  t      =  cosd(lon0) .* u + sind(lon0) .* v;
  uEast  = -sind(lon0) .* u + cosd(lon0) .* v;
  wUp    =  cosd(lat0) .* t + sind(lat0) .* w;
  vNorth = -sind(lat0) .* t + cosd(lat0) .* w;
end
