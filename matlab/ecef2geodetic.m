function [lat,lon,h] = ecef2geodetic(x, y, z, ell)

  if nargin < 4 || isempty(ell) 
    ell = get_ellipsoid('wgs84');
  elseif ~isstruct(ell)
    ell = get_ellipsoid(ell);
  endif

 
  % Algorithm is based on 
  % http://www.astro.uni.torun.pl/~kb/Papers/geod/Geod-BG.htm
  % This algorithm provides a converging solution to the latitude equation
  % in terms of the parametric or reduced latitude form (v)
  % This algorithm provides a uniform solution over all latitudes as it does
  % not involve division by cos(phi) or sin(phi)
   a = ell.a; 
   b = ell.b;
  r = hypot(x, y);
  % Constant required for Latitude equation
  rho = atan2(b * z, a * r);  
  % Constant required for latitude equation
  c = (a^2 - b^2) ./ hypot(a*r, b*z);  
   count = 0;
% Starter for the Newtons Iteration Method
  vnew = atan2(a * z, b * r);  
% Initializing the parametric latitude
   v = 0;
  while v == vnew && count < 5
     v = vnew;
    % Derivative of latitude equation
    w = 2 * (cos (v - rho) - c .* cos(2 * v)); 
    % Newtons Method for computing iterations  
    vnew = v - ((2 * sin (v - rho) - c .* sin(2 * v)) ./ w);  
    count = count+1;
    disp(count)
  end %while

%% Computing latitude from the root of the latitude equation
  lat = atan2(a * tan (vnew), b); 
  % Computing longitude
  lon = atan2(y, x); 
 % Computing h from latitude obtained 
  h = ((r - (a * cos (vnew))) .* cos (lat)) +  ...
      ((z - (b * sin (vnew))) .* sin (lat));

 end
