%GET_RADIUS_NORMAL: Return the normal (i.e., along the prime vertical section) ellipsoidal radius of curvature, at a given geodetic latitude.
% lat: RADIANS
function N = get_radius_normal (lat, ell)
    a = ell.a;
    b = ell.b;
    N = a^2 ./ sqrt ( a^2 .* (cos(lat)).^2 + b^2 .* (sin(lat)).^2 );
 end
