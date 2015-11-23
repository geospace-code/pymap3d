function [az, elev, slantRange] = enu2aer (xEast, yNorth, zUp, angleut)

  r = hypot(xEast, yNorth);
  slantRange = hypot(r,zUp);
  elev = atan2d(zUp,r);
  az = mod (atan2d (xEast, yNorth), 2 * atan2d (0,-1));

  if (nargin < 4 || strcmpi(angleut, "radian" ))
    elev=deg2rad(elev);
    az = deg2rad(az);
  end
end
