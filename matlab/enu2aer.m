function [az, elev, slantRange] = enu2aer(e, n, u, angleut)

  
  

  r = hypot(e, n);
  slantRange = hypot(r,u);
  % radians
  elev = atan2(u,r);
  az = mod(atan2(e, n), 2 * atan2(0,-1));

  if nargin < 4 || strcmpi(angleut(1),'d')
    elev = rad2deg(elev);
    az = rad2deg(az);
  end
  
end
