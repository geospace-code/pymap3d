"""Transforms involving NED North East Down
They are typically just reordered and sign-changed versions of ENU functions,
included for convenience and readability
"""

from __future__ import annotations

from ._typing import FloatLike

from .ecef import ecef2enu, ecef2enuv, ecef2geodetic, enu2ecef
from .ellipsoid import Ellipsoid
from .enu import aer2enu, enu2aer, geodetic2enu
from .frames import _ecef2ned_rotation, _ned2ecef_rotation

__all__ = [
    "aer2ned",
    "ned2aer",
    "ned2geodetic",
    "ned2ecef",
    "ecef2ned",
    "geodetic2ned",
    "ecef2nedv",
    "ecef2ned_matrix",
    "ned2ecef_matrix",
]


def aer2ned(az: FloatLike, elev: FloatLike, slantRange: FloatLike, deg: bool = True) -> tuple:
    """
    converts azimuth, elevation, range to target from observer to North, East, Down

    Parameters
    -----------

    az : array-like float
        azimuth
    elev : array-like float
        elevation
    slantRange : array-like float
        slant range [meters]
    deg : bool, optional
        degrees input/output  (False: radians in/out)

    Results
    -------
    n : array-like float
        North NED coordinate (meters)
    e : array-like float
        East NED coordinate (meters)
    d : array-like float
        Down NED coordinate (meters)
    """
    e, n, u = aer2enu(az, elev, slantRange, deg=deg)

    return n, e, -u


def ned2aer(n: FloatLike, e: FloatLike, d: FloatLike, deg: bool = True) -> tuple:
    """
    converts North, East, Down to azimuth, elevation, range

    Parameters
    ----------

    n : array-like float
        North NED coordinate (meters)
    e : array-like float
        East NED coordinate (meters)
    d : array-like float
        Down NED coordinate (meters)
    deg : bool, optional
        degrees input/output  (False: radians in/out)

    Results
    -------

    az : array-like float
        azimuth
    elev : array-like float
        elevation
    slantRange : array-like float
         slant range [meters]
    """

    return enu2aer(e, n, -d, deg=deg)


def ned2geodetic(
    n: FloatLike,
    e: FloatLike,
    d: FloatLike,
    lat0: FloatLike,
    lon0: FloatLike,
    h0: FloatLike,
    ell: Ellipsoid | None = None,
    deg: bool = True,
) -> tuple:
    """
    Converts North, East, Down to target latitude, longitude, altitude

    Parameters
    ----------

    n : array-like float
        North NED coordinate (meters)
    e : array-like float
        East NED coordinate (meters)
    d : array-like float
        Down NED coordinate (meters)
    lat0 : array-like float
        Observer geodetic latitude
    lon0 : array-like float
        Observer geodetic longitude
    h0 : array-like float
        Observer altitude above geodetic ellipsoid (meters)
    ell : Ellipsoid, optional
        reference ellipsoid
    deg : bool, optional
        degrees input/output  (False: radians in/out)

    Results
    -------

    lat : array-like float
        target geodetic latitude
    lon : array-like float
        target geodetic longitude
    h : array-like float
        target altitude above geodetic ellipsoid (meters)
    """
    x, y, z = enu2ecef(e, n, -d, lat0, lon0, h0, ell, deg=deg)

    return ecef2geodetic(x, y, z, ell, deg=deg)


def ned2ecef(
    n: FloatLike,
    e: FloatLike,
    d: FloatLike,
    lat0: FloatLike,
    lon0: FloatLike,
    h0: FloatLike,
    ell: Ellipsoid | None = None,
    deg: bool = True,
) -> tuple:
    """
    North, East, Down to target ECEF coordinates

    Parameters
    ----------

    n : array-like float
        North NED coordinate (meters)
    e : array-like float
        East NED coordinate (meters)
    d : array-like float
        Down NED coordinate (meters)
    lat0 : array-like float
        Observer geodetic latitude
    lon0 : array-like float
        Observer geodetic longitude
    h0 : array-like float
        Observer altitude above geodetic ellipsoid (meters)
    ell : Ellipsoid, optional
        reference ellipsoid
    deg : bool, optional
        degrees input/output  (False: radians in/out)

    Results
    -------

    x : array-like float
        ECEF x coordinate (meters)
    y : array-like float
        ECEF y coordinate (meters)
    z : array-like float
        ECEF z coordinate (meters)
    """
    return enu2ecef(e, n, -d, lat0, lon0, h0, ell, deg=deg)


def ecef2ned(
    x: FloatLike,
    y: FloatLike,
    z: FloatLike,
    lat0: FloatLike,
    lon0: FloatLike,
    h0: FloatLike,
    ell: Ellipsoid | None = None,
    deg: bool = True,
) -> tuple:
    """
    Convert ECEF x,y,z to North, East, Down

    Parameters
    ----------

    x : array-like float
        ECEF x coordinate (meters)
    y : array-like float
        ECEF y coordinate (meters)
    z : array-like float
        ECEF z coordinate (meters)
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

    n : array-like float
        North NED coordinate (meters)
    e : array-like float
        East NED coordinate (meters)
    d : array-like float
        Down NED coordinate (meters)
    """
    e, n, u = ecef2enu(x, y, z, lat0, lon0, h0, ell, deg=deg)

    return n, e, -u


def geodetic2ned(
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
    convert latitude, longitude, altitude of target to North, East, Down from observer

    Parameters
    ----------

    lat : array-like float
        target geodetic latitude
    lon : array-like float
        target geodetic longitude
    h : array-like float
        target altitude above geodetic ellipsoid (meters)
    lat0 : array-like float
        Observer geodetic latitude
    lon0 : array-like float
        Observer geodetic longitude
    h0 : array-like float
        Observer altitude above geodetic ellipsoid (meters)
    ell : Ellipsoid, optional
        reference ellipsoid
    deg : bool, optional
        degrees input/output  (False: radians in/out)


    Results
    -------

    n : array-like float
        North NED coordinate (meters)
    e : array-like float
        East NED coordinate (meters)
    d : array-like float
        Down NED coordinate (meters)
    """
    e, n, u = geodetic2enu(lat, lon, h, lat0, lon0, h0, ell, deg=deg)

    return n, e, -u


def ecef2nedv(
    x: FloatLike, y: FloatLike, z: FloatLike, lat0: FloatLike, lon0: FloatLike, deg: bool = True
) -> tuple:
    """
    for VECTOR between two points

    Parameters
    ----------
    x : array-like float
        ECEF x coordinate (meters)
    y : array-like float
        ECEF y coordinate (meters)
    z : array-like float
        ECEF z coordinate (meters)
    lat0 : array-like float
        Observer geodetic latitude
    lon0 : array-like float
        Observer geodetic longitude
    deg : bool, optional
          degrees input/output  (False: radians in/out)

    Results
    -------

    (Vector)

    n : array-like float
        North NED coordinate (meters)
    e : array-like float
        East NED coordinate (meters)
    d : array-like float
        Down NED coordinate (meters)
    """
    e, n, u = ecef2enuv(x, y, z, lat0, lon0, deg=deg)

    return n, e, -u


def ecef2ned_matrix(lat0, lon0, deg: bool = True):
    """
    Rotation matrix that maps ECEF vectors into NED coordinates.
    """

    return _ecef2ned_rotation(lat0, lon0, deg=deg)


def ned2ecef_matrix(lat0, lon0, deg: bool = True):
    """
    Rotation matrix that maps NED vectors into ECEF coordinates.
    """

    return _ned2ecef_rotation(lat0, lon0, deg=deg)
