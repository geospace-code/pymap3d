"""transforms involving ENU East North Up"""

from __future__ import annotations

from math import tau

try:
    from numpy import asarray
except ImportError:
    pass

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


def enu2aer(e, n, u, deg: bool = True) -> tuple:
    """
    ENU to Azimuth, Elevation, Range

    Parameters
    ----------

    e : float
        ENU East coordinate (meters)
    n : float
        ENU North coordinate (meters)
    u : float
        ENU Up coordinate (meters)
    deg : bool, optional
        degrees input/output  (False: radians in/out)

    Results
    -------

    azimuth : float
        azimuth to rarget
    elevation : float
        elevation to target
    srange : float
        slant range [meters]
    """

    # 1 millimeter precision for singularity stability

    try:
        e[abs(e) < 1e-3] = 0.0
        n[abs(n) < 1e-3] = 0.0
        u[abs(u) < 1e-3] = 0.0
    except TypeError:
        if abs(e) < 1e-3:
            e = 0.0
        if abs(n) < 1e-3:
            n = 0.0
        if abs(u) < 1e-3:
            u = 0.0

    r = hypot(e, n)
    slantRange = hypot(r, u)
    elev = atan2(u, r)
    az = atan2(e, n) % tau

    if deg:
        az = degrees(az)
        elev = degrees(elev)

    return az, elev, slantRange


def aer2enu(az, el, srange, deg: bool = True) -> tuple:
    """
    Azimuth, Elevation, Slant range to target to East, North, Up

    Parameters
    ----------
    az : float
            azimuth clockwise from north (degrees)
    el : float
        elevation angle above horizon, neglecting aberrations (degrees)
    srange : float
        slant range [meters]
    deg : bool, optional
        degrees input/output  (False: radians in/out)

    Returns
    --------
    e : float
        East ENU coordinate (meters)
    n : float
        North ENU coordinate (meters)
    u : float
        Up ENU coordinate (meters)
    """
    if deg:
        el = radians(el)
        az = radians(az)

    try:
        if (asarray(srange) < 0).any():
            raise ValueError("Slant range  [0, Infinity)")
    except NameError:
        if srange < 0:
            raise ValueError("Slant range  [0, Infinity)")

    r = srange * cos(el)

    return r * sin(az), r * cos(az), srange * sin(el)


def enu2geodetic(
    e,
    n,
    u,
    lat0,
    lon0,
    h0,
    ell: Ellipsoid | None = None,
    deg: bool = True,
) -> tuple:
    """
    East, North, Up to target to geodetic coordinates

    Parameters
    ----------
    e : float
        East ENU coordinate (meters)
    n : float
        North ENU coordinate (meters)
    u : float
        Up ENU coordinate (meters)
    lat0 : float
           Observer geodetic latitude
    lon0 : float
           Observer geodetic longitude
    h0 : float
         observer altitude above geodetic ellipsoid (meters)
    ell : Ellipsoid, optional
          reference ellipsoid
    deg : bool, optional
          degrees input/output  (False: radians in/out)


    Results
    -------
    lat : float
          geodetic latitude
    lon : float
          geodetic longitude
    alt : float
          altitude above ellipsoid  (meters)
    """

    x, y, z = enu2ecef(e, n, u, lat0, lon0, h0, ell, deg=deg)

    return ecef2geodetic(x, y, z, ell, deg=deg)


def geodetic2enu(
    lat,
    lon,
    h,
    lat0,
    lon0,
    h0,
    ell: Ellipsoid | None = None,
    deg: bool = True,
) -> tuple:
    """
    Parameters
    ----------
    lat : float
          target geodetic latitude
    lon : float
          target geodetic longitude
    h : float
          target altitude above ellipsoid  (meters)
    lat0 : float
           Observer geodetic latitude
    lon0 : float
           Observer geodetic longitude
    h0 : float
         observer altitude above geodetic ellipsoid (meters)
    ell : Ellipsoid, optional
          reference ellipsoid
    deg : bool, optional
          degrees input/output  (False: radians in/out)


    Results
    -------
    e : float
        East ENU
    n : float
        North ENU
    u : float
        Up ENU
    """
    x1, y1, z1 = geodetic2ecef(lat, lon, h, ell, deg=deg)
    x2, y2, z2 = geodetic2ecef(lat0, lon0, h0, ell, deg=deg)

    return uvw2enu(x1 - x2, y1 - y2, z1 - z2, lat0, lon0, deg=deg)


def enu2ecefv(e, n, u, lat0, lon0, deg: bool = True) -> tuple:
    """
    VECTOR from observer to target ENU => ECEF

    Parameters
    ----------
    e
        target e ENU coordinate
    n
        target n ENU coordinate
    u
        target u ENU coordinate
    lat0
           Observer geodetic latitude
    lon0
           Observer geodetic longitude
    deg : bool, optional
          degrees input/output  (False: radians in/out)

    Returns
    -------
    x
        target x ECEF coordinate
    y
        target y ECEF coordinate
    z
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
