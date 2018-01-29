% Copyright (c) 2014-2018 Michael Hirsch, Ph.D.
% Copyright (c) 2013, Felipe Geremia Nievinski
%
% Redistribution and use in source and binary forms, with or without modification, are permitted provided that the following conditions are met:
% 1. Redistributions of source code must retain the above copyright notice, this list of conditions and the following disclaimer.
% 2. Redistributions in binary form must reproduce the above copyright notice, this list of conditions and the following disclaimer in the documentation and/or other materials provided with the distribution.
% THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.


function [e, n, Up] = ecef2enuv (u, v, w, lat0, lon0, angleUnit)
% function [e, n, Up] = ecef2enuv (u, v, w, lat0, lon0, angleUnit)
%
% Inputs
% ------
% u,v,w: meters
% lat0,lon0: geodetic latitude and longitude (degrees)
% angleUnit: string for angular system. Default 'd' degrees
%
% Outputs
% -------
% e,n,Up:  East, North, Up vector
%
  if nargin<6 || isempty(angleUnit) || strcmpi(angleUnit(1), 'd')
    lat0 = deg2rad(lat0);
    lon0 = deg2rad(lon0);
  end

  t  =  cos(lon0) .* u + sin(lon0) .* v;
  e  = -sin(lon0) .* u + cos(lon0) .* v;
  Up =  cos(lat0) .* t + sin(lat0) .* w;
  n  = -sin(lat0) .* t + cos(lat0) .* w;
end
