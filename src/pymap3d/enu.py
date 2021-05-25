""" transforms involving ENU East North Up """
from __future__ import annotations
import typing

from math import tau

try:
    from numpy import asarray, radians, sin, cos, hypot, arctan2 as atan2, degrees, ndarray
except ImportError:
    from math import radians, sin, cos, hypot, atan2, degrees  # type: ignore

    ndarray = typing.Any  # type: ignore

from .ecef import geodetic2ecef, ecef2geodetic, enu2ecef, uvw2enu
from .ellipsoid import Ellipsoid

__all__ = ["enu2aer", "aer2enu", "enu2geodetic", "geodetic2enu"]


def enu2aer(
    e: ndarray, n: ndarray, u: ndarray, deg: bool = True
) -> tuple[ndarray, ndarray, ndarray]:
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

    # 1 millimeter precision for singularity

    try:
        e[abs(e) < 1e-3] = 0.0
        n[abs(n) < 1e-3] = 0.0
        u[abs(u) < 1e-3] = 0.0
    except TypeError:
        if abs(e) < 1e-3:
            e = 0.0  # type: ignore
        if abs(n) < 1e-3:
            n = 0.0  # type: ignore
        if abs(u) < 1e-3:
            u = 0.0  # type: ignore

    r = hypot(e, n)
    slantRange = hypot(r, u)
    elev = atan2(u, r)
    az = atan2(e, n) % tau

    if deg:
        az = degrees(az)
        elev = degrees(elev)

    return az, elev, slantRange


def aer2enu(
    az: ndarray, el: ndarray, srange: float | ndarray, deg: bool = True
) -> tuple[ndarray, ndarray, ndarray]:
    """
    Azimuth, Elevation, Slant range to target to East, North, Up

    Parameters
    ----------
    azimuth : float
            azimuth clockwise from north (degrees)
    elevation : float
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
    e: ndarray,
    n: ndarray,
    u: ndarray,
    lat0: ndarray,
    lon0: ndarray,
    h0: ndarray,
    ell: Ellipsoid = None,
    deg: bool = True,
) -> tuple[ndarray, ndarray, ndarray]:
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
    lat: ndarray,
    lon: ndarray,
    h: ndarray,
    lat0: ndarray,
    lon0: ndarray,
    h0: ndarray,
    ell: Ellipsoid = None,
    deg: bool = True,
) -> tuple[ndarray, ndarray, ndarray]:
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
