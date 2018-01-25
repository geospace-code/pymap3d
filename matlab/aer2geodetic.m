function [lat1, lon1, h1] = aer2geodetic (az, el, slantRange, lat0, lon0, ...
                            h0, ell, angleut)
                            
  if nargin<8
    angleut='d';
  end

  %% Convert AER2ECEF
  [x, y, z] = aer2ecef (az, el, slantRange, lat0, lon0, h0, ell, angleut);
  %% Convert ECEF2GEODETIC
  [lat1, lon1, h1] = ecef2geodetic (x, y, z, ell, angleut);

end
