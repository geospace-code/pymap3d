"""transforms involving ENU East North Up"""

from __future__ import annotations

from math import tau

from ._typing import FloatLike

from .ecef import ecef2geodetic, enu2ecef, geodetic2ecef, uvw2enu
from .ellipsoid import Ellipsoid
from .frames import _ecef2enu_rotation, _enu2ecef_rotation
from .mathfun import atan2, cos, degrees, hypot, radians, sin

__all__ = [
    "enu2aer",
    "aer2enu",
    "enu2geodetic",
    "geodetic2enu",
    "enu2ecefv",
    "ecef2enu_matrix",
    "enu2ecef_matrix",
]


def enu2aer(e: FloatLike, n: FloatLike, u: FloatLike, deg: bool = True) -> tuple:
    """
    ENU to Azimuth, Elevation, Range

    Parameters
    ----------

    e : array-like float
        ENU East coordinate (meters)
    n : array-like float
        ENU North coordinate (meters)
    u : array-like float
        ENU Up coordinate (meters)
    deg : bool, optional
        degrees input/output  (False: radians in/out)

    Results
    -------

    azimuth : array-like float
        azimuth to rarget
    elevation : array-like float
        elevation to target
    srange : array-like float
        slant range [meters]
    """

    # 1 millimeter precision for singularity stability

    def _zero_trim_numpy(x, y, z) -> tuple:
        x[abs(x) < 1e-3] = 0.0
        y[abs(y) < 1e-3] = 0.0
        z[abs(z) < 1e-3] = 0.0
        return x, y, z

    def _zero_trim_scalar(*values: float, threshold: float = 1e-3) -> tuple[float, ...]:
        """Trim small values (|x| < threshold) to zero."""
        return tuple(0.0 if abs(v) < threshold else float(v) for v in values)

    try:
        e, n, u = _zero_trim_numpy(e, n, u)
    except TypeError:
        e, n, u = _zero_trim_scalar(float(e), float(n), float(u))

    r = hypot(e, n)
    slantRange = hypot(r, u)
    elev = atan2(u, r)
    az = atan2(e, n) % tau

    if deg:
        az = degrees(az)
        elev = degrees(elev)

    return az, elev, slantRange


def aer2enu(az: FloatLike, el: FloatLike, srange: FloatLike, deg: bool = True) -> tuple:
    """
    Azimuth, Elevation, Slant range to target to East, North, Up

    Parameters
    ----------
    az : array-like float
        azimuth clockwise from north (degrees)
    el : array-like float
        elevation angle above horizon, neglecting aberrations (degrees)
    srange : array-like float
        slant range [meters]. expected to be non-negative.
    deg : bool, optional
        degrees input/output  (False: radians in/out)

    Returns
    --------
    e : array-like float
        East ENU coordinate (meters)
    n : array-like float
        North ENU coordinate (meters)
    u : array-like float
        Up ENU coordinate (meters)
    """

    if deg:
        el = radians(el)
        az = radians(az)

    r = srange * cos(el)

    return r * sin(az), r * cos(az), srange * sin(el)


def enu2geodetic(
    e: FloatLike,
    n: FloatLike,
    u: FloatLike,
    lat0: FloatLike,
    lon0: FloatLike,
    h0: FloatLike,
    ell: Ellipsoid | None = None,
    deg: bool = True,
) -> tuple:
    """
    East, North, Up to target to geodetic coordinates

    Parameters
    ----------
    e : array-like float
        East ENU coordinate (meters)
    n : array-like float
        North ENU coordinate (meters)
    u : array-like float
        Up ENU coordinate (meters)
    lat0 : array-like float
        Observer geodetic latitude
    lon0 : array-like float
        Observer geodetic longitude
    h0 : array-like float
        observer altitude above geodetic ellipsoid (meters)
    ell : Ellipsoid, optional
        reference ellipsoid
    deg : bool, optional
        degrees input/output  (False: radians in/out)


    Results
    -------
    lat : array-like float
          geodetic latitude
    lon : array-like float
          geodetic longitude
    alt : array-like float
          altitude above ellipsoid  (meters)
    """

    x, y, z = enu2ecef(e, n, u, lat0, lon0, h0, ell, deg=deg)

    return ecef2geodetic(x, y, z, ell, deg=deg)


def geodetic2enu(
    lat: FloatLike,
    lon: FloatLike,
    h: FloatLike,
    lat0: FloatLike,
    lon0: FloatLike,
    h0: FloatLike,
    ell: Ellipsoid | None = None,
    deg: bool = True,
) -> tuple:
    """
    Parameters
    ----------
    lat : array-like float
          target geodetic latitude
    lon : array-like float
          target geodetic longitude
    h : array-like float
          target altitude above ellipsoid  (meters)
    lat0 : array-like float
           Observer geodetic latitude
    lon0 : array-like float
           Observer geodetic longitude
    h0 : array-like float
         observer altitude above geodetic ellipsoid (meters)
    ell : Ellipsoid, optional
          reference ellipsoid
    deg : bool, optional
          degrees input/output  (False: radians in/out)


    Results
    -------
    e : array-like float
        East ENU
    n : array-like float
        North ENU
    u : array-like float
        Up ENU
    """
    x1, y1, z1 = geodetic2ecef(lat, lon, h, ell, deg=deg)
    x2, y2, z2 = geodetic2ecef(lat0, lon0, h0, ell, deg=deg)

    return uvw2enu(x1 - x2, y1 - y2, z1 - z2, lat0, lon0, deg=deg)


def enu2ecefv(
    e: FloatLike, n: FloatLike, u: FloatLike, lat0: FloatLike, lon0: FloatLike, deg: bool = True
) -> tuple:
    """
    VECTOR from observer to target ENU => ECEF

    Parameters
    ----------
    e : array-like float
        target e ENU coordinate
    n : array-like float
        target n ENU coordinate
    u : array-like float
        target u ENU coordinate
    lat0 : array-like float
        Observer geodetic latitude
    lon0 : array-like float
        Observer geodetic longitude
    deg : bool, optional
        degrees input/output  (False: radians in/out)

    Returns
    -------
    x : array-like float
        target x ECEF coordinate
    y : array-like float
        target y ECEF coordinate
    z : array-like float
        target z ECEF coordinate

    """
    if deg:
        lat0 = radians(lat0)
        lon0 = radians(lon0)

    x = -sin(lon0) * e - sin(lat0) * cos(lon0) * n + cos(lat0) * cos(lon0) * u
    y = cos(lon0) * e - sin(lat0) * sin(lon0) * n + cos(lat0) * sin(lon0) * u
    z = cos(lat0) * n + sin(lat0) * u

    return x, y, z


def ecef2enu_matrix(lat0, lon0, deg: bool = True):
    """
    Rotation matrix that maps ECEF vectors into ENU coordinates.
    """

    return _ecef2enu_rotation(lat0, lon0, deg=deg)


def enu2ecef_matrix(lat0, lon0, deg: bool = True):
    """
    Rotation matrix that maps ENU vectors into ECEF coordinates.
    """

    return _enu2ecef_rotation(lat0, lon0, deg=deg)
