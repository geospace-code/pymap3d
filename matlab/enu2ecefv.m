% Copyright (c) 2014-2018 Michael Hirsch, Ph.D.
% Copyright (c) 2013, Felipe Geremia Nievinski
%
% Redistribution and use in source and binary forms, with or without modification, are permitted provided that the following conditions are met:
% 1. Redistributions of source code must retain the above copyright notice, this list of conditions and the following disclaimer.
% 2. Redistributions in binary form must reproduce the above copyright notice, this list of conditions and the following disclaimer in the documentation and/or other materials provided with the distribution.
% THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.


function [u, v, w] = enu2ecefv(uEast, vNorth, wUp, lat0, lon0, angleut)

  if nargin < 6 
      angleut = 'd';
  end
      
  if strcmpi(angleut(1),'r')
    lat0 = rad2deg(lat0);
    lon0 = rad2deg(lon0);
  end

  t = cos(lat0) .* wUp - sin(lat0) .* vNorth;
  w = sin(lat0) .* wUp + cos(lat0) .* vNorth;

  u = cos(lon0) .* t - sin(lon0) .* uEast;
  v = sin(lon0) .* t + cos(lon0) .* uEast;
  
end
