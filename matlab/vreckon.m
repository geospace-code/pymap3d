function [lat2,lon2,a21] = vreckon(lat1,lon1,s,a12)
% RECKON - Using the WGS-84 Earth ellipsoid, travel a given distance along
%          a given azimuth starting at a given initial point, and return
%          the endpoint within a few millimeters of accuracy, using
%          Vincenty's algorithm.
%
% USAGE:
% [lat2,lon2] = vreckon(lat1, lon1, s, a12)
%
% VARIABLES:
% lat1 = inital latitude (degrees)
% lon1 = initial longitude (degrees)
% s    = distance (meters)
% a12  = intial azimuth (degrees)
% lat2, lon2 = second point (degrees)
% a21  = reverse azimuth (degrees), at final point facing back toward the
%        intial point
%
% Original algorithm source:
% T. Vincenty, "Direct and Inverse Solutions of Geodesics on the Ellipsoid
% with Application of Nested Equations", Survey Review, vol. 23, no. 176,
% April 1975, pp 88-93.
% Available at: http://www.ngs.noaa.gov/PUBS_LIB/inverse.pdf
%
% Notes: 
% (1) The Vincenty reckoning algorithm was transcribed verbatim into
%     JavaScript by Chris Veness. It was modified and translated to Matlab
%     by Michael Kleder. Mr. Veness's website is:
%     http://www.movable-type.co.uk/scripts/latlong-vincenty-direct.html
% (2) Error correcting code, polar error corrections, WGS84 ellipsoid
%     parameters, testing, and comments by Michael Kleder.
% (3) By convention, when starting at a pole, the longitude of the initial
%     point (otherwise meaningless) determines the longitude line along
%     which to traverse, and hence the longitude of the final point.
% (4) The convention noted in (3) above creates a discrepancy with VDIST
%     when the the intial or final point is at a pole. In the VDIST
%     function, when traversing from a pole, the azimuth is  0 when
%     heading away from the south pole and 180 when heading away from the
%     north pole. In contrast, this VRECKON function uses the azimuth as
%     noted in (3) above when traversing away form a pole.
% (5) In testing, where the traversal subtends no more than 178 degrees,
%     this function correctly inverts the VDIST function to within 0.2
%     millimeters of distance, 5e-10 degrees of forward azimuth,
%     and 5e-10 degrees of reverse azimuth. Precision reduces as test
%     points approach antipodal because the precision of VDIST is reduced
%     for nearly antipodal points. (A warning is given by VDIST.)
% (6) Tested but no warranty. Use at your own risk.
% (7) Ver 1.0, Michael Kleder, November 2007

% Input check:
if abs(lat1)>90
    error('Input latitude must be between -90 and 90 degrees, inclusive.')
end
a = 6378137; % semimajor axis
b = 6356752.31424518; % semiminor axis
f = 1/298.257223563; % flattening coefficient
lat1   = lat1 * .1745329251994329577e-1; % intial latitude in radians
lon1   = lon1 * .1745329251994329577e-1; % intial longitude in radians
% correct for errors at exact poles by adjusting 0.6 millimeters:
kidx = abs(pi/2-abs(lat1)) < 1e-10;
if any(kidx);
    lat1(kidx) = sign(lat1(kidx))*(pi/2-(1e-10));
end
alpha1 = a12 * .1745329251994329577e-1; % inital azimuth in radians
sinAlpha1 = sin(alpha1);
cosAlpha1 = cos(alpha1);
tanU1 = (1-f) * tan(lat1);
cosU1 = 1 / sqrt(1 + tanU1*tanU1);
sinU1 = tanU1*cosU1;
sigma1 = atan2(tanU1, cosAlpha1);
sinAlpha = cosU1 * sinAlpha1;
cosSqAlpha = 1 - sinAlpha*sinAlpha;
uSq = cosSqAlpha * (a*a - b*b) / (b*b);
A = 1 + uSq/16384*(4096+uSq*(-768+uSq*(320-175*uSq)));
B = uSq/1024 * (256+uSq*(-128+uSq*(74-47*uSq)));
sigma = s / (b*A);
sigmaP = 2*pi;
while (abs(sigma-sigmaP) > 1e-12)
    cos2SigmaM = cos(2*sigma1 + sigma);
    sinSigma = sin(sigma);
    cosSigma = cos(sigma);
    deltaSigma = B*sinSigma*(cos2SigmaM+B/4*(cosSigma*(-1+...
        2*cos2SigmaM*cos2SigmaM)-...
        B/6*cos2SigmaM*(-3+4*sinSigma*sinSigma)*(-3+...
        4*cos2SigmaM*cos2SigmaM)));
    sigmaP = sigma;
    sigma = s / (b*A) + deltaSigma;
end
tmp = sinU1*sinSigma - cosU1*cosSigma*cosAlpha1;
lat2 = atan2(sinU1*cosSigma + cosU1*sinSigma*cosAlpha1,...
    (1-f)*sqrt(sinAlpha*sinAlpha + tmp*tmp));
lambda = atan2(sinSigma*sinAlpha1, cosU1*cosSigma - ...
    sinU1*sinSigma*cosAlpha1);
C = f/16*cosSqAlpha*(4+f*(4-3*cosSqAlpha));
L = lambda - (1-C) * f * sinAlpha * (sigma + C*sinSigma*(cos2SigmaM+...
    C*cosSigma*(-1+2*cos2SigmaM*cos2SigmaM)));
lon2 = lon1 + L;
% output degrees
lat2 = lat2 * 57.295779513082322865;
lon2 = lon2 * 57.295779513082322865;
lon2 = mod(lon2,360); % follow [0,360] convention
if nargout > 2
    a21 = atan2(sinAlpha, -tmp); 
    a21  = 180 + a21  * 57.295779513082322865; % note direction reversal
    a21=mod(a21,360);
end
return