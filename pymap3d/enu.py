""" transforms involving ENU East North Up """
from typing import Tuple
from numpy import radians, sin, cos, hypot, arctan2, degrees
import numpy as np
try:
    from math import tau
except ImportError:
    tau = 2 * np.pi

from .ecef import geodetic2ecef, ecef2geodetic, enu2ecef, uvw2enu


__all__ = ['enu2aer', 'aer2enu', 'enu2geodetic', 'geodetic2enu']


def enu2aer(e: np.ndarray, n: np.ndarray, u: np.ndarray, deg: bool = True) -> Tuple[float, float, float]:
    """
    ENU to Azimuth, Elevation, Range

    Parameters
    ----------

    e : float or np.ndarray of float
        ENU East coordinate (meters)
    n : float or np.ndarray of float
        ENU North coordinate (meters)
    u : float or np.ndarray of float
        ENU Up coordinate (meters)
    deg : bool, optional
        degrees input/output  (False: radians in/out)

    Results
    -------

    azimuth : float or np.ndarray of float
        azimuth to rarget
    elevation : float or np.ndarray of float
        elevation to target
    srange : float or np.ndarray of float
        slant range [meters]
    """
    # 1 millimeter precision for singularity

    e = np.asarray(e)
    n = np.asarray(n)
    u = np.asarray(u)

    with np.errstate(invalid='ignore'):
        e[abs(e) < 1e-3] = 0.
        n[abs(n) < 1e-3] = 0.
        u[abs(u) < 1e-3] = 0.

        r = hypot(e, n)
        slantRange = hypot(r, u)
        elev = arctan2(u, r)
        az = arctan2(e, n) % tau

    if deg:
        az = degrees(az)
        elev = degrees(elev)

    return az, elev, slantRange


def aer2enu(az: float, el: float, srange: float, deg: bool = True) -> Tuple[float, float, float]:
    """
    Azimuth, Elevation, Slant range to target to East, north, Up

    Parameters
    ----------
    azimuth : float or np.ndarray of float
            azimuth clockwise from north (degrees)
    elevation : float or np.ndarray of float
        elevation angle above horizon, neglecting aberattions (degrees)
    srange : float or np.ndarray of float
        slant range [meters]
    deg : bool, optional
        degrees input/output  (False: radians in/out)

    Returns
    --------
    e : float or np.ndarray of float
        East ENU coordinate (meters)
    n : float or np.ndarray of float
        North ENU coordinate (meters)
    u : float or np.ndarray of float
        Up ENU coordinate (meters)
    """
    if deg:
        el = radians(el)
        az = radians(az)

    with np.errstate(invalid='ignore'):
        if (np.asarray(srange) < 0).any():
            raise ValueError('Slant range  [0, Infinity)')

    r = srange * cos(el)

    return r * sin(az), r * cos(az), srange * sin(el)


def enu2geodetic(e: float, n: float, u: float,
                 lat0: float, lon0: float, h0: float,
                 ell=None, deg: bool = True) -> Tuple[float, float, float]:
    """
    East, North, Up to target to geodetic coordinates

    Parameters
    ----------
    e : float or np.ndarray of float
        East ENU coordinate (meters)
    n : float or np.ndarray of float
        North ENU coordinate (meters)
    u : float or np.ndarray of float
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
    lat : float or np.ndarray of float
          geodetic latitude
    lon : float or np.ndarray of float
          geodetic longitude
    alt : float or np.ndarray of float
          altitude above ellipsoid  (meters)
    """

    x, y, z = enu2ecef(e, n, u, lat0, lon0, h0, ell, deg=deg)

    return ecef2geodetic(x, y, z, ell, deg=deg)


def geodetic2enu(lat: float, lon: float, h: float,
                 lat0: float, lon0: float, h0: float,
                 ell=None, deg: bool = True) -> Tuple[float, float, float]:
    """
    Parameters
    ----------
    lat : float or np.ndarray of float
          target geodetic latitude
    lon : float or np.ndarray of float
          target geodetic longitude
    h : float or np.ndarray of float
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
    e : float or np.ndarray of float
        East ENU
    n : float or np.ndarray of float
        North ENU
    u : float or np.ndarray of float
        Up ENU
    """
    x1, y1, z1 = geodetic2ecef(lat, lon, h, ell, deg=deg)
    x2, y2, z2 = geodetic2ecef(lat0, lon0, h0, ell, deg=deg)
    dx = x1 - x2
    dy = y1 - y2
    dz = z1 - z2

    return uvw2enu(dx, dy, dz, lat0, lon0, deg=deg)
