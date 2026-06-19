from __future__ import annotations

from .mathfun import degrees, radians, sin, cos, atan2, asin
from .ecef import geodetic2ecef, ecef2geodetic
from .ellipsoid import Ellipsoid


def geodetic2nvector(lat, lon, deg: bool = True) -> tuple:
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


def nvector2geodetic(n1, n2, n3, deg=True) -> tuple:
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


def ecef2nvector(x, y, z, ell: Ellipsoid | None = None, deg: bool = True):
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


def nvector2ecef(n1, n2, n3, alt=0, ell: Ellipsoid | None = None, deg: bool = True):
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
