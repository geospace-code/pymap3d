function varargout = vdist(lat1,lon1,lat2,lon2)
% VDIST - Using the WGS-84 Earth ellipsoid, compute the distance between
%         two points within a few millimeters of accuracy, compute forward
%         azimuth, and compute backward azimuth, all using a vectorized
%         version of Vincenty's algorithm.
%
% s = vdist(lat1,lon1,lat2,lon2)
% [s,a12] = vdist(lat1,lon1,lat2,lon2)
% [s,a12,a21] = vdist(lat1,lon1,lat2,lon2)
%
% s = distance in meters (inputs may be scalars, vectors, or matrices)
% a12 = azimuth in degrees from first point to second point (forward)
% a21 = azimuth in degrees from second point to first point (backward)
%       (Azimuths are in degrees clockwise from north.)
% lat1 = GEODETIC latitude of first point (degrees)
% lon1 = longitude of first point (degrees)
% lat2, lon2 = second point (degrees)
%
%  Original algorithm source:
%  T. Vincenty, "Direct and Inverse Solutions of Geodesics on the Ellipsoid
%  with Application of Nested Equations", Survey Review, vol. 23, no. 176,
%  April 1975, pp 88-93.
%  Available at: http://www.ngs.noaa.gov/PUBS_LIB/inverse.pdf
%
% Notes: (1) lat1,lon1,lat2,lon2 can be any (identical) size/shape. Outputs
%            will have the same size and shape.
%        (2) Error correcting code, convergence failure traps, antipodal
%            corrections, polar error corrections, WGS84 ellipsoid
%            parameters, testing, and comments: Michael Kleder, 2004.
%        (3) Azimuth implementation (including quadrant abiguity
%            resolution) and code vectorization, Michael Kleder, Sep 2005.
%        (4) Vectorization is convergence sensitive; that is, quantities
%            which have already converged to within tolerance are not
%            recomputed during subsequent iterations (while other
%            quantities are still converging).
%        (5) Vincenty describes his distance algorithm as precise to within
%            0.01 millimeters, subject to the ellipsoidal model.
%        (6) For distance calculations, essentially antipodal points are
%            treated as exactly antipodal, potentially reducing accuracy
%            slightly.
%        (7) Distance failures for points exactly at the poles are
%            eliminated by moving the points by 0.6 millimeters.
%        (8) The Vincenty distance algorithm was transcribed verbatim by
%            Peter Cederholm, August 12, 2003. It was modified and
%            translated to English by Michael Kleder.
%            Mr. Cederholm's website is http://www.plan.aau.dk/~pce/
%        (9) Distances agree with the Mapping Toolbox, version 2.2 (R14SP3)
%            with a max relative difference of about 5e-9, except when the
%            two points are nearly antipodal, and except when one point is
%            near the equator and the two longitudes are nearly 180 degrees
%            apart. This function (vdist) is more accurate in such cases.
%            For example, note this difference (as of this writing):
%            >>vdist(0.2,305,15,125)
%            18322827.0131551
%            >>distance(0.2,305,15,125,[6378137 0.08181919])
%            0
%       (10) Azimuths FROM the north pole (either forward starting at the
%            north pole or backward when ending at the north pole) are set
%            to 180 degrees by convention. Azimuths FROM the south pole are
%            set to 0 degrees by convention.
%       (11) Azimuths agree with the Mapping Toolbox, version 2.2 (R14SP3)
%            to within about a hundred-thousandth of a degree, except when
%            traversing to or from a pole, where the convention for this
%            function is described in (10), and except in the cases noted
%            above in (9).
%       (12) No warranties; use at your own risk.

% reshape inputs
keepsize = size(lat1);
lat1=lat1(:);
lon1=lon1(:);
lat2=lat2(:);
lon2=lon2(:);
% Input check:
if any(abs(lat1)>90 | abs(lat2)>90)
    error('Input latitudes must be between -90 and 90 degrees, inclusive.')
end
% Supply WGS84 earth ellipsoid axis lengths in meters:
a = 6378137; % definitionally
b = 6356752.31424518; % computed from WGS84 earth flattening coefficient
% preserve true input latitudes:
lat1tr = lat1;
lat2tr = lat2;
% convert inputs in degrees to radians:
lat1 = lat1 * 0.0174532925199433;
lon1 = lon1 * 0.0174532925199433;
lat2 = lat2 * 0.0174532925199433;
lon2 = lon2 * 0.0174532925199433;
% correct for errors at exact poles by adjusting 0.6 millimeters:
kidx = abs(pi/2-abs(lat1)) < 1e-10;
if any(kidx);
    lat1(kidx) = sign(lat1(kidx))*(pi/2-(1e-10));
end
kidx = abs(pi/2-abs(lat2)) < 1e-10;
if any(kidx)
    lat2(kidx) = sign(lat2(kidx))*(pi/2-(1e-10));
end
f = (a-b)/a;
U1 = atan((1-f)*tan(lat1));
U2 = atan((1-f)*tan(lat2));
lon1 = mod(lon1,2*pi);
lon2 = mod(lon2,2*pi);
L = abs(lon2-lon1);
kidx = L > pi;
if any(kidx)
    L(kidx) = 2*pi - L(kidx);
end
lambda = L;
lambdaold = 0*lat1;
itercount = 0;
notdone = logical(1+0*lat1);
alpha = 0*lat1;
sigma = 0*lat1;
cos2sigmam = 0*lat1;
C = 0*lat1;
warninggiven = logical(0);
while any(notdone)  % force at least one execution
    %disp(['lambda(21752) = ' num2str(lambda(21752),20)]);
    itercount = itercount+1;
    if itercount > 50
        if ~warninggiven
            warning(['Essentially antipodal points encountered. ' ...
                'Precision may be reduced slightly.']);
        end
        lambda(notdone) = pi;
        break
    end
    lambdaold(notdone) = lambda(notdone);
    sinsigma(notdone) = sqrt((cos(U2(notdone)).*sin(lambda(notdone)))...
        .^2+(cos(U1(notdone)).*sin(U2(notdone))-sin(U1(notdone)).*...
        cos(U2(notdone)).*cos(lambda(notdone))).^2);
    cossigma(notdone) = sin(U1(notdone)).*sin(U2(notdone))+...
        cos(U1(notdone)).*cos(U2(notdone)).*cos(lambda(notdone));
    % eliminate rare imaginary portions at limit of numerical precision:
    sinsigma(notdone)=real(sinsigma(notdone));
    cossigma(notdone)=real(cossigma(notdone));
    sigma(notdone) = atan2(sinsigma(notdone),cossigma(notdone));
    alpha(notdone) = asin(cos(U1(notdone)).*cos(U2(notdone)).*...
        sin(lambda(notdone))./sin(sigma(notdone)));
    cos2sigmam(notdone) = cos(sigma(notdone))-2*sin(U1(notdone)).*...
        sin(U2(notdone))./cos(alpha(notdone)).^2;
    C(notdone) = f/16*cos(alpha(notdone)).^2.*(4+f*(4-3*...
        cos(alpha(notdone)).^2));
    lambda(notdone) = L(notdone)+(1-C(notdone)).*f.*sin(alpha(notdone))...
        .*(sigma(notdone)+C(notdone).*sin(sigma(notdone)).*...
        (cos2sigmam(notdone)+C(notdone).*cos(sigma(notdone)).*...
        (-1+2.*cos2sigmam(notdone).^2)));
    %disp(['then, lambda(21752) = ' num2str(lambda(21752),20)]);
    % correct for convergence failure in the case of essentially antipodal
    % points
    if any(lambda(notdone) > pi)
        warning(['Essentially antipodal points encountered. ' ...
            'Precision may be reduced slightly.']);
        warninggiven = logical(1);
        lambdaold(lambda>pi) = pi;
        lambda(lambda>pi) = pi;
    end
    notdone = abs(lambda-lambdaold) > 1e-12;
end
u2 = cos(alpha).^2.*(a^2-b^2)/b^2;
A = 1+u2./16384.*(4096+u2.*(-768+u2.*(320-175.*u2)));
B = u2./1024.*(256+u2.*(-128+u2.*(74-47.*u2)));
deltasigma = B.*sin(sigma).*(cos2sigmam+B./4.*(cos(sigma).*(-1+2.*...
    cos2sigmam.^2)-B./6.*cos2sigmam.*(-3+4.*sin(sigma).^2).*(-3+4*...
    cos2sigmam.^2)));
varargout{1} = reshape(b.*A.*(sigma-deltasigma),keepsize);
if nargout > 1
    % From point #1 to point #2
    % correct sign of lambda for azimuth calcs:
    lambda = abs(lambda);
    kidx=sign(sin(lon2-lon1)) .* sign(sin(lambda)) < 0;
    lambda(kidx) = -lambda(kidx);
    numer = cos(U2).*sin(lambda);
    denom = cos(U1).*sin(U2)-sin(U1).*cos(U2).*cos(lambda);
    a12 = atan2(numer,denom);
    kidx = a12<0;
    a12(kidx)=a12(kidx)+2*pi;
    % from poles:
    a12(lat1tr <= -90) = 0;
    a12(lat1tr >= 90 ) = pi;
    varargout{2} = reshape(a12 * 57.2957795130823,keepsize); % to degrees
end
if nargout > 2
    a21=NaN*lat1;
    % From point #2 to point #1
    % correct sign of lambda for azimuth calcs:
    lambda = abs(lambda);
    kidx=sign(sin(lon1-lon2)) .* sign(sin(lambda)) < 0;
    lambda(kidx)=-lambda(kidx);
    numer = cos(U1).*sin(lambda);
    denom = sin(U1).*cos(U2)-cos(U1).*sin(U2).*cos(lambda);
    a21 = atan2(numer,denom);
    kidx=a21<0;
    a21(kidx)= a21(kidx)+2*pi;
    % backwards from poles:
    a21(lat2tr >= 90) = pi;
    a21(lat2tr <= -90) = 0;
    varargout{3} = reshape(a21 * 57.2957795130823,keepsize); % to degrees
end
return

