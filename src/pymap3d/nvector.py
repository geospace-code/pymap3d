"""
N-vector based geodetic operations.

Reference
---------
Gade, K. (2010). A Nonsingular Horizontal Position Representation,
The Journal of Navigation, Volume 63, Issue 03, pp 395-417, July 2010.
(www.navlab.net/Publications/A_Nonsingular_Horizontal_Position_Representation.pdf)
"""

from __future__ import annotations

from ._typing import FloatLike
from .mathfun import acos, asin, atan2, cos, degrees, radians, sin, sqrt
from .ecef import geodetic2ecef, ecef2geodetic
from .ellipsoid import Ellipsoid


def geodetic2nvector(lat: FloatLike, lon: FloatLike, deg: bool = True) -> tuple:
    """
    Convert geodetic coordinates (latitude, longitude) to an n-vector.

    Parameters:
        lat : array-like float
            Geodetic latitude(s).
        lon : array-like float
            Geodetic longitude(s).
        ell : str or tuple, optional
            Reference ellipsoid (default is None, which uses WGS84).
        deg : bool, optional
            If True (default), inputs are in degrees. If False, use radians.

    Returns:
        n1, n2, n3 : ndarray
            Components of the n-vector in the Earth-Centered Earth-Fixed (ECEF) coordinate system.
    """

    if deg:
        lat, lon = radians(lat), radians(lon)

    sin_lat, cos_lat = sin(lat), cos(lat)
    sin_lon, cos_lon = sin(lon), cos(lon)

    n1 = cos_lat * cos_lon
    n2 = cos_lat * sin_lon
    n3 = sin_lat

    return n1, n2, n3


def nvector2geodetic(n1: FloatLike, n2: FloatLike, n3: FloatLike, deg: bool = True) -> tuple:
    """
    Convert an n-vector back to geodetic coordinates (latitude, longitude).

    Parameters:
        n1, n2, n3 : array-like float
            Components of the n-vector in the Earth-Centered Earth-Fixed (ECEF) coordinate system.
        ell : str or tuple, optional
            Reference ellipsoid (default is None, which uses WGS84).
        deg : bool, optional
            If True (default), returns latitude and longitude in degrees. If False, in radians.

    Returns:
        lat, lon : ndarray
            Geodetic latitude(s) and longitude(s).
    """

    # Compute latitude and longitude from n-vector
    lat = asin(n3)
    lon = atan2(n2, n1)

    if deg:
        lat, lon = degrees(lat), degrees(lon)

    return lat, lon


def ecef2nvector(
    x: FloatLike, y: FloatLike, z: FloatLike, ell: Ellipsoid | None = None, deg: bool = True
):
    """
    Convert ECEF coordinates to an n-vector.

    Parameters:
        x, y, z : array-like float
            ECEF coordinates in meters.
        ell : str or tuple, optional
            Reference ellipsoid (default is None, which uses WGS84).
        deg : bool, optional
            If True (default), geodetic2nvector() inputs are in degrees. If False, in radians.

    Returns:
        n1, n2, n3 : ndarray
            Components of the n-vector in the Earth-Centered Earth-Fixed (ECEF) coordinate system.
    """

    lat, lon, _ = ecef2geodetic(x, y, z, ell=ell, deg=deg)
    return geodetic2nvector(lat, lon, deg=deg)


def nvector2ecef(
    n1: FloatLike,
    n2: FloatLike,
    n3: FloatLike,
    alt: FloatLike = 0,
    ell: Ellipsoid | None = None,
    deg: bool = True,
):
    """
    Convert an n-vector to ECEF coordinates.

    Parameters:
        n1, n2, n3 : array-like float
            Components of the n-vector in the Earth-Centered Earth-Fixed (ECEF) coordinate system.
        alt : array-like float, optional
            Altitude in meters (default is 0).
        ell : str or tuple, optional
            Reference ellipsoid (default is None, which uses WGS84).
        deg : bool, optional
            If True (default), nvector2geodetic() outputs are in degrees. If False, in radians.

    Returns:
        x, y, z : ndarray
            ECEF coordinates in meters.
    """
    lat, lon = nvector2geodetic(n1, n2, n3, deg=deg)
    return geodetic2ecef(lat, lon, alt, ell=ell, deg=deg)


# --- Private helpers for 3-component vector math ---


def _dot3(a1, a2, a3, b1, b2, b3):
    """Dot product of two 3-vectors given as components."""
    return a1 * b1 + a2 * b2 + a3 * b3


def _cross3(a1, a2, a3, b1, b2, b3):
    """Cross product of two 3-vectors given as components."""
    return (a2 * b3 - a3 * b2, a3 * b1 - a1 * b3, a1 * b2 - a2 * b1)


def _norm3(a1, a2, a3):
    """Normalize a 3-vector to unit length."""
    r = sqrt(a1 * a1 + a2 * a2 + a3 * a3)
    return (a1 / r, a2 / r, a3 / r)


# --- N-vector spherical geometry operations ---


def _default_radius():
    """Mean Earth radius (equal-volume sphere) in meters."""
    from .rsphere import eqavol

    return eqavol()


def nvector_distance(n1_1, n1_2, n1_3, n2_1, n2_2, n2_3, radius: float | None = None) -> float:
    """
    Great-circle distance between two n-vectors.

    Uses atan2 formulation for numerical stability at all angular separations.

    Parameters
    ----------
    n1_1, n1_2, n1_3 : float or array-like
        Components of first n-vector.
    n2_1, n2_2, n2_3 : float or array-like
        Components of second n-vector.
    radius : float, optional
        Sphere radius in meters (default: WGS84 equal-volume sphere, ~6371 km).

    Returns
    -------
    dist : float
        Great-circle distance in meters.
    """
    if radius is None:
        radius = _default_radius()

    cx, cy, cz = _cross3(n1_1, n1_2, n1_3, n2_1, n2_2, n2_3)
    cross_mag = sqrt(cx * cx + cy * cy + cz * cz)
    dot = _dot3(n1_1, n1_2, n1_3, n2_1, n2_2, n2_3)

    return radius * atan2(cross_mag, dot)


def nvector_interpolate(n1_1, n1_2, n1_3, n2_1, n2_2, n2_3, fraction: float) -> tuple:
    """
    Spherical linear interpolation (SLERP) between two n-vectors.

    Parameters
    ----------
    n1_1, n1_2, n1_3 : float
        Components of first n-vector.
    n2_1, n2_2, n2_3 : float
        Components of second n-vector.
    fraction : float
        Interpolation fraction. 0.0 returns n1, 1.0 returns n2.

    Returns
    -------
    ni_1, ni_2, ni_3 : float
        Components of interpolated n-vector.
    """
    dot = _dot3(n1_1, n1_2, n1_3, n2_1, n2_2, n2_3)
    # clamp to [-1, 1] for safety
    dot = max(-1.0, min(1.0, dot))

    omega = acos(dot)

    # For very small angles, fall back to normalized linear interpolation
    if abs(omega) < 1e-10:
        i1 = (1.0 - fraction) * n1_1 + fraction * n2_1
        i2 = (1.0 - fraction) * n1_2 + fraction * n2_2
        i3 = (1.0 - fraction) * n1_3 + fraction * n2_3
        return _norm3(i1, i2, i3)

    s = sin(omega)
    a = sin((1.0 - fraction) * omega) / s
    b = sin(fraction * omega) / s

    i1 = a * n1_1 + b * n2_1
    i2 = a * n1_2 + b * n2_2
    i3 = a * n1_3 + b * n2_3

    return _norm3(i1, i2, i3)


def nvector_mean(n1s, n2s, n3s) -> tuple:
    """
    Mean geographic position of multiple n-vectors.

    Computes the normalized centroid of the input n-vectors on the unit sphere.

    Parameters
    ----------
    n1s, n2s, n3s : array-like
        Sequences of n-vector components.

    Returns
    -------
    n1, n2, n3 : float
        Components of the mean n-vector (normalized to unit length).
    """
    s1 = sum(n1s)
    s2 = sum(n2s)
    s3 = sum(n3s)

    return _norm3(s1, s2, s3)


def nvector_cross_track_distance(
    n1_1, n1_2, n1_3,
    n2_1, n2_2, n2_3,
    np_1, np_2, np_3,
    radius: float | None = None,
) -> float:
    """
    Cross-track distance from a point to a great circle path.

    Parameters
    ----------
    n1_1, n1_2, n1_3 : float
        Components of n-vector for first point on the path.
    n2_1, n2_2, n2_3 : float
        Components of n-vector for second point on the path.
    np_1, np_2, np_3 : float
        Components of n-vector for the query point.
    radius : float, optional
        Sphere radius in meters (default: WGS84 equal-volume sphere).

    Returns
    -------
    dist : float
        Signed cross-track distance in meters. Positive means the point
        is to the left of the path (from n1 toward n2).
    """
    if radius is None:
        radius = _default_radius()

    # Great circle pole: normal to the plane containing n1 and n2
    cx, cy, cz = _cross3(n1_1, n1_2, n1_3, n2_1, n2_2, n2_3)
    cx, cy, cz = _norm3(cx, cy, cz)

    # Signed angular distance from the great circle
    return radius * asin(_dot3(cx, cy, cz, np_1, np_2, np_3))


def nvector_intersection(
    n1_1, n1_2, n1_3,
    n2_1, n2_2, n2_3,
    n3_1, n3_2, n3_3,
    n4_1, n4_2, n4_3,
) -> tuple:
    """
    Intersection of two great circle paths.

    Path 1 passes through (n1, n2), path 2 passes through (n3, n4).
    Returns one of the two antipodal intersection points, chosen as the
    one closest to the centroid of all four input points.

    Parameters
    ----------
    n1_1, n1_2, n1_3 : float
        Components of first n-vector on path 1.
    n2_1, n2_2, n2_3 : float
        Components of second n-vector on path 1.
    n3_1, n3_2, n3_3 : float
        Components of first n-vector on path 2.
    n4_1, n4_2, n4_3 : float
        Components of second n-vector on path 2.

    Returns
    -------
    ni_1, ni_2, ni_3 : float
        Components of the intersection n-vector.
    """
    # Normals to each great circle plane
    c1_1, c1_2, c1_3 = _cross3(n1_1, n1_2, n1_3, n2_1, n2_2, n2_3)
    c2_1, c2_2, c2_3 = _cross3(n3_1, n3_2, n3_3, n4_1, n4_2, n4_3)

    # Intersection is perpendicular to both normals
    i1, i2, i3 = _cross3(c1_1, c1_2, c1_3, c2_1, c2_2, c2_3)
    i1, i2, i3 = _norm3(i1, i2, i3)

    # Pick the solution closest to the centroid of all four points
    m1 = n1_1 + n2_1 + n3_1 + n4_1
    m2 = n1_2 + n2_2 + n3_2 + n4_2
    m3 = n1_3 + n2_3 + n3_3 + n4_3

    if _dot3(i1, i2, i3, m1, m2, m3) < 0:
        i1, i2, i3 = -i1, -i2, -i3

    return i1, i2, i3
