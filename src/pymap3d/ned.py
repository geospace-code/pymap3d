""" Transforms involving NED North East Down """

import typing

from .enu import geodetic2enu, aer2enu, enu2aer
from .ecef import ecef2geodetic, ecef2enuv, ecef2enu, enu2ecef
from .ellipsoid import Ellipsoid

try:
    from numpy.typing import ArrayLike
except ImportError:
    ArrayLike = typing.Any


def aer2ned(
    az: ArrayLike, elev: ArrayLike, slantRange: ArrayLike, deg: bool = True
) -> typing.Tuple[ArrayLike, ArrayLike, ArrayLike]:
    """
    converts azimuth, elevation, range to target from observer to North, East, Down

    Parameters
    -----------

    az : ArrayLike
         azimuth
    elev : ArrayLike
         elevation
    slantRange : ArrayLike
         slant range [meters]
    deg : bool, optional
          degrees input/output  (False: radians in/out)

    Results
    -------
    n : ArrayLike
        North NED coordinate (meters)
    e : ArrayLike
        East NED coordinate (meters)
    d : ArrayLike
        Down NED coordinate (meters)
    """
    e, n, u = aer2enu(az, elev, slantRange, deg=deg)

    return n, e, -u


def ned2aer(n: ArrayLike, e: ArrayLike, d: ArrayLike, deg: bool = True) -> typing.Tuple[ArrayLike, ArrayLike, ArrayLike]:
    """
    converts North, East, Down to azimuth, elevation, range

    Parameters
    ----------

    n : ArrayLike
        North NED coordinate (meters)
    e : ArrayLike
        East NED coordinate (meters)
    d : ArrayLike
        Down NED coordinate (meters)
    deg : bool, optional
        degrees input/output  (False: radians in/out)

    Results
    -------

    az : ArrayLike
         azimuth
    elev : ArrayLike
         elevation
    slantRange : ArrayLike
         slant range [meters]
    """
    return enu2aer(e, n, -d, deg=deg)


def ned2geodetic(
    n: ArrayLike,
    e: ArrayLike,
    d: ArrayLike,
    lat0: ArrayLike,
    lon0: ArrayLike,
    h0: ArrayLike,
    ell: Ellipsoid = None,
    deg: bool = True,
) -> typing.Tuple[ArrayLike, ArrayLike, ArrayLike]:
    """
    Converts North, East, Down to target latitude, longitude, altitude

    Parameters
    ----------

    n : ArrayLike
        North NED coordinate (meters)
    e : ArrayLike
        East NED coordinate (meters)
    d : ArrayLike
        Down NED coordinate (meters)
    lat0 : ArrayLike
        Observer geodetic latitude
    lon0 : ArrayLike
        Observer geodetic longitude
    h0 : ArrayLike
         observer altitude above geodetic ellipsoid (meters)
    ell : Ellipsoid, optional
          reference ellipsoid
    deg : bool, optional
          degrees input/output  (False: radians in/out)

    Results
    -------

    lat : ArrayLike
        target geodetic latitude
    lon : ArrayLike
        target geodetic longitude
    h : ArrayLike
        target altitude above geodetic ellipsoid (meters)

    """
    x, y, z = enu2ecef(e, n, -d, lat0, lon0, h0, ell, deg=deg)

    return ecef2geodetic(x, y, z, ell, deg=deg)


def ned2ecef(
    n: ArrayLike,
    e: ArrayLike,
    d: ArrayLike,
    lat0: ArrayLike,
    lon0: ArrayLike,
    h0: ArrayLike,
    ell: Ellipsoid = None,
    deg: bool = True,
) -> typing.Tuple[ArrayLike, ArrayLike, ArrayLike]:
    """
    North, East, Down to target ECEF coordinates

    Parameters
    ----------

    n : ArrayLike
        North NED coordinate (meters)
    e : ArrayLike
        East NED coordinate (meters)
    d : ArrayLike
        Down NED coordinate (meters)
    lat0 : ArrayLike
        Observer geodetic latitude
    lon0 : ArrayLike
        Observer geodetic longitude
    h0 : ArrayLike
         observer altitude above geodetic ellipsoid (meters)
    ell : Ellipsoid, optional
          reference ellipsoid
    deg : bool, optional
          degrees input/output  (False: radians in/out)

    Results
    -------

    x : ArrayLike
        ECEF x coordinate (meters)
    y : ArrayLike
        ECEF y coordinate (meters)
    z : ArrayLike
        ECEF z coordinate (meters)
    """
    return enu2ecef(e, n, -d, lat0, lon0, h0, ell, deg=deg)


def ecef2ned(
    x: ArrayLike,
    y: ArrayLike,
    z: ArrayLike,
    lat0: ArrayLike,
    lon0: ArrayLike,
    h0: ArrayLike,
    ell: Ellipsoid = None,
    deg: bool = True,
) -> typing.Tuple[ArrayLike, ArrayLike, ArrayLike]:
    """
    Convert ECEF x,y,z to North, East, Down

    Parameters
    ----------

    x : ArrayLike
        ECEF x coordinate (meters)
    y : ArrayLike
        ECEF y coordinate (meters)
    z : ArrayLike
        ECEF z coordinate (meters)
    lat0 : ArrayLike
        Observer geodetic latitude
    lon0 : ArrayLike
        Observer geodetic longitude
    h0 : ArrayLike
         observer altitude above geodetic ellipsoid (meters)
    ell : Ellipsoid, optional
          reference ellipsoid
    deg : bool, optional
          degrees input/output  (False: radians in/out)

    Results
    -------

    n : ArrayLike
        North NED coordinate (meters)
    e : ArrayLike
        East NED coordinate (meters)
    d : ArrayLike
        Down NED coordinate (meters)

    """
    e, n, u = ecef2enu(x, y, z, lat0, lon0, h0, ell, deg=deg)

    return n, e, -u


def geodetic2ned(
    lat: ArrayLike,
    lon: ArrayLike,
    h: ArrayLike,
    lat0: ArrayLike,
    lon0: ArrayLike,
    h0: ArrayLike,
    ell: Ellipsoid = None,
    deg: bool = True,
) -> typing.Tuple[ArrayLike, ArrayLike, ArrayLike]:
    """
    convert latitude, longitude, altitude of target to North, East, Down from observer

    Parameters
    ----------

    lat : ArrayLike
        target geodetic latitude
    lon : ArrayLike
        target geodetic longitude
    h : ArrayLike
        target altitude above geodetic ellipsoid (meters)
    lat0 : ArrayLike
        Observer geodetic latitude
    lon0 : ArrayLike
        Observer geodetic longitude
    h0 : ArrayLike
         observer altitude above geodetic ellipsoid (meters)
    ell : Ellipsoid, optional
          reference ellipsoid
    deg : bool, optional
          degrees input/output  (False: radians in/out)


    Results
    -------

    n : ArrayLike
        North NED coordinate (meters)
    e : ArrayLike
        East NED coordinate (meters)
    d : ArrayLike
        Down NED coordinate (meters)
    """
    e, n, u = geodetic2enu(lat, lon, h, lat0, lon0, h0, ell, deg=deg)

    return n, e, -u


def ecef2nedv(
    x: ArrayLike, y: ArrayLike, z: ArrayLike, lat0: ArrayLike, lon0: ArrayLike, deg: bool = True
) -> typing.Tuple[ArrayLike, ArrayLike, ArrayLike]:
    """
    for VECTOR between two points

    Parameters
    ----------
    x : ArrayLike
        ECEF x coordinate (meters)
    y : ArrayLike
        ECEF y coordinate (meters)
    z : ArrayLike
        ECEF z coordinate (meters)
    lat0 : ArrayLike
        Observer geodetic latitude
    lon0 : ArrayLike
        Observer geodetic longitude
    deg : bool, optional
          degrees input/output  (False: radians in/out)

    Results
    -------

    (Vector)

    n : ArrayLike
        North NED coordinate (meters)
    e : ArrayLike
        East NED coordinate (meters)
    d : ArrayLike
        Down NED coordinate (meters)
    """
    e, n, u = ecef2enuv(x, y, z, lat0, lon0, deg=deg)

    return n, e, -u
