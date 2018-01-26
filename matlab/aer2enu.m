function [e, n, u] = aer2enu (az, el, slantRange, angleut)

  if nargin==3 || strcmpi(angleut(1),'d') 
    az = deg2rad(az);
    el = deg2rad(el);
  end    

%% Calculation of AER2ENU
   u = slantRange .* sin(el);
   r = slantRange .* cos(el);
   e = r .* sin(az);
   n = r .* cos(az);

end
