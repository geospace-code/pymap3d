function [lat1, lon1, h1] = aer2geodetic (az, el, slantRange, lat0, lon0, ...
                            h0, ell, angleut)
     
  if nargin<7 || isempty(ell), ell=get_ellipsoid(); end
  if nargin<8, angleut='d';  end

  [x, y, z] = aer2ecef(az, el, slantRange, lat0, lon0, h0, ell, angleut);

  [lat1, lon1, h1] = ecef2geodetic(x, y, z, ell, angleut);

end
