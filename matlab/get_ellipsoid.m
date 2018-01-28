% Copyright (c) 2014-2018 Michael Hirsch, Ph.D.
% Copyright (c) 2013, Felipe Geremia Nievinski
%
% Redistribution and use in source and binary forms, with or without modification, are permitted provided that the following conditions are met:
% 1. Redistributions of source code must retain the above copyright notice, this list of conditions and the following disclaimer.
% 2. Redistributions in binary form must reproduce the above copyright notice, this list of conditions and the following disclaimer in the documentation and/or other materials provided with the distribution.
% THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

function ell = get_ellipsoid(model) 
% function ell = get_ellipsoid(model) 
%
% inputs
% ------
% model: string of model name e.g. 'wgs84'
%
%
% outputs
% -------
% ell: struct containing ellipsoid parameters
%
%
if ~nargin
  model='wgs84';
end

switch model
case 'wgs84'
 % WGS-84 ellipsoid parameters.
 %   http://earth-info.nga.mil/GandG/tr8350_2.html
 %   ftp://164.214.2.65/pub/gig/tr8350.2/wgs84fin.pdf
 ell.a = 6378137.0;                                % semi-major axis
 ell.f = 1.0 / 298.2572235630;                     % flattening
 ell.b = ell.a * (1 - ell.f);                      % semi-minor axis
 %ell.e = sqrt ( (ell.a^2 - ell.b^2) / (ell.a^2));  % first eccentricity
case 'grs80'
% GRS-80 ellipsoid parameters
% <http://itrf.ensg.ign.fr/faq.php?type=answer> (accessed 2018-01-22)
 ell.a = 6378137.0;                                % semi-major axis
 ell.f = 1/298.257222100882711243;                 % flattening
 ell.b = ell.a * (1 - ell.f);                      % semi-minor axis
 %ell.e = sqrt ( (ell.a^2 - ell.b^2) / (ell.a^2));  % first eccentricity 
otherwise, error([name,' not yet implemented'])
end

end % function
