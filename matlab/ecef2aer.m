% Copyright (c) 2014-2018 Michael Hirsch, Ph.D.
% Copyright (c) 2013, Felipe Geremia Nievinski
%
% Redistribution and use in source and binary forms, with or without modification, are permitted provided that the following conditions are met:
% 1. Redistributions of source code must retain the above copyright notice, this list of conditions and the following disclaimer.
% 2. Redistributions in binary form must reproduce the above copyright notice, this list of conditions and the following disclaimer in the documentation and/or other materials provided with the distribution.
% THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.


function [az, el, slantRange] = ecef2aer (x, y, z, lat0, lon0, alt0, spheroid, angleUnit)
% function [az, el, slantRange] = ecef2aer (x, y, z, lat0, lon0, alt0, spheroid, angleUnit)
%
% Inputs
% ------
% x,y,z: Earth Centered Earth Fixed (ECEF) coordinates of test point (meters)
%
% Outputs
% -------
% az, el, slantrange: look angles and distance to point under test (degrees, degrees, meters)
% az: azimuth clockwise from local north
% el: elevation angle above local horizon

  if nargin < 7 || isempty (ell), spheroid = wgs84Ellipsoid();  end
  if nargin<8, angleUnit='d'; end

  [e, n, u] = ecef2enu(x, y, z, lat0, lon0, alt0, spheroid, angleUnit);
  [az,el,slantRange] = enu2aer(e, n, u, angleUnit);
  
end
