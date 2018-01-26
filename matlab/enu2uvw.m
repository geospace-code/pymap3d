function [u,v,w] = enu2uvw(e,n,u,lat0,lon0,angleut)
  
    if nargin<6 || strcmpi(angleut(1),'d')
        lat0 = deg2rad(lat0)
        lon0 = deg2rad(lon0)
    end
    
    t = cos(lat0) * u - sin(lat0) * n
    w = sin(lat0) * u + cos(lat0) * n

    u = cos(lon0) * t - sin(lon0) * e
    v = sin(lon0) * t + cos(lon0) * e
  
end % function