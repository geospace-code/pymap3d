""" Transforms involving NED North East Down """

import typing

from .enu import geodetic2enu, aer2enu, enu2aer
from .ecef import ecef2geodetic, ecef2enuv, ecef2enu, enu2ecef
from .ellipsoid import Ellipsoid

if typing.TYPE_CHECKING:
    from numpy import ndarray


def aer2ned(
    az: "ndarray", elev: "ndarray", slantRange: "ndarray", deg: bool = True
) -> typing.Tuple["ndarray", "ndarray", "ndarray"]:
    """
    converts azimuth, elevation, range to target from observer to North, East, Down

    Parameters
    -----------

    az : "ndarray"
         azimuth
    elev : "ndarray"
         elevation
    slantRange : "ndarray"
         slant range [meters]
    deg : bool, optional
          degrees input/output  (False: radians in/out)

    Results
    -------
    n : "ndarray"
        North NED coordinate (meters)
    e : "ndarray"
        East NED coordinate (meters)
    d : "ndarray"
        Down NED coordinate (meters)
    """
    e, n, u = aer2enu(az, elev, slantRange, deg=deg)

    return n, e, -u


def ned2aer(n: "ndarray", e: "ndarray", d: "ndarray", deg: bool = True) -> typing.Tuple["ndarray", "ndarray", "ndarray"]:
    """
    converts North, East, Down to azimuth, elevation, range

    Parameters
    ----------

    n : "ndarray"
        North NED coordinate (meters)
    e : "ndarray"
        East NED coordinate (meters)
    d : "ndarray"
        Down NED coordinate (meters)
    deg : bool, optional
        degrees input/output  (False: radians in/out)

    Results
    -------

    az : "ndarray"
         azimuth
    elev : "ndarray"
         elevation
    slantRange : "ndarray"
         slant range [meters]
    """
    return enu2aer(e, n, -d, deg=deg)


def ned2geodetic(
    n: "ndarray",
    e: "ndarray",
    d: "ndarray",
    lat0: "ndarray",
    lon0: "ndarray",
    h0: "ndarray",
    ell: Ellipsoid = None,
    deg: bool = True,
) -> typing.Tuple["ndarray", "ndarray", "ndarray"]:
    """
    Converts North, East, Down to target latitude, longitude, altitude

    Parameters
    ----------

    n : "ndarray"
        North NED coordinate (meters)
    e : "ndarray"
        East NED coordinate (meters)
    d : "ndarray"
        Down NED coordinate (meters)
    lat0 : "ndarray"
        Observer geodetic latitude
    lon0 : "ndarray"
        Observer geodetic longitude
    h0 : "ndarray"
         observer altitude above geodetic ellipsoid (meters)
    ell : Ellipsoid, optional
          reference ellipsoid
    deg : bool, optional
          degrees input/output  (False: radians in/out)

    Results
    -------

    lat : "ndarray"
        target geodetic latitude
    lon : "ndarray"
        target geodetic longitude
    h : "ndarray"
        target altitude above geodetic ellipsoid (meters)

    """
    x, y, z = enu2ecef(e, n, -d, lat0, lon0, h0, ell, deg=deg)

    return ecef2geodetic(x, y, z, ell, deg=deg)


def ned2ecef(
    n: "ndarray",
    e: "ndarray",
    d: "ndarray",
    lat0: "ndarray",
    lon0: "ndarray",
    h0: "ndarray",
    ell: Ellipsoid = None,
    deg: bool = True,
) -> typing.Tuple["ndarray", "ndarray", "ndarray"]:
    """
    North, East, Down to target ECEF coordinates

    Parameters
    ----------

    n : "ndarray"
        North NED coordinate (meters)
    e : "ndarray"
        East NED coordinate (meters)
    d : "ndarray"
        Down NED coordinate (meters)
    lat0 : "ndarray"
        Observer geodetic latitude
    lon0 : "ndarray"
        Observer geodetic longitude
    h0 : "ndarray"
         observer altitude above geodetic ellipsoid (meters)
    ell : Ellipsoid, optional
          reference ellipsoid
    deg : bool, optional
          degrees input/output  (False: radians in/out)

    Results
    -------

    x : "ndarray"
        ECEF x coordinate (meters)
    y : "ndarray"
        ECEF y coordinate (meters)
    z : "ndarray"
        ECEF z coordinate (meters)
    """
    return enu2ecef(e, n, -d, lat0, lon0, h0, ell, deg=deg)


def ecef2ned(
    x: "ndarray",
    y: "ndarray",
    z: "ndarray",
    lat0: "ndarray",
    lon0: "ndarray",
    h0: "ndarray",
    ell: Ellipsoid = None,
    deg: bool = True,
) -> typing.Tuple["ndarray", "ndarray", "ndarray"]:
    """
    Convert ECEF x,y,z to North, East, Down

    Parameters
    ----------

    x : "ndarray"
        ECEF x coordinate (meters)
    y : "ndarray"
        ECEF y coordinate (meters)
    z : "ndarray"
        ECEF z coordinate (meters)
    lat0 : "ndarray"
        Observer geodetic latitude
    lon0 : "ndarray"
        Observer geodetic longitude
    h0 : "ndarray"
         observer altitude above geodetic ellipsoid (meters)
    ell : Ellipsoid, optional
          reference ellipsoid
    deg : bool, optional
          degrees input/output  (False: radians in/out)

    Results
    -------

    n : "ndarray"
        North NED coordinate (meters)
    e : "ndarray"
        East NED coordinate (meters)
    d : "ndarray"
        Down NED coordinate (meters)

    """
    e, n, u = ecef2enu(x, y, z, lat0, lon0, h0, ell, deg=deg)

    return n, e, -u


def geodetic2ned(
    lat: "ndarray",
    lon: "ndarray",
    h: "ndarray",
    lat0: "ndarray",
    lon0: "ndarray",
    h0: "ndarray",
    ell: Ellipsoid = None,
    deg: bool = True,
) -> typing.Tuple["ndarray", "ndarray", "ndarray"]:
    """
    convert latitude, longitude, altitude of target to North, East, Down from observer

    Parameters
    ----------

    lat : "ndarray"
        target geodetic latitude
    lon : "ndarray"
        target geodetic longitude
    h : "ndarray"
        target altitude above geodetic ellipsoid (meters)
    lat0 : "ndarray"
        Observer geodetic latitude
    lon0 : "ndarray"
        Observer geodetic longitude
    h0 : "ndarray"
         observer altitude above geodetic ellipsoid (meters)
    ell : Ellipsoid, optional
          reference ellipsoid
    deg : bool, optional
          degrees input/output  (False: radians in/out)


    Results
    -------

    n : "ndarray"
        North NED coordinate (meters)
    e : "ndarray"
        East NED coordinate (meters)
    d : "ndarray"
        Down NED coordinate (meters)
    """
    e, n, u = geodetic2enu(lat, lon, h, lat0, lon0, h0, ell, deg=deg)

    return n, e, -u


def ecef2nedv(
    x: "ndarray", y: "ndarray", z: "ndarray", lat0: "ndarray", lon0: "ndarray", deg: bool = True
) -> typing.Tuple["ndarray", "ndarray", "ndarray"]:
    """
    for VECTOR between two points

    Parameters
    ----------
    x : "ndarray"
        ECEF x coordinate (meters)
    y : "ndarray"
        ECEF y coordinate (meters)
    z : "ndarray"
        ECEF z coordinate (meters)
    lat0 : "ndarray"
        Observer geodetic latitude
    lon0 : "ndarray"
        Observer geodetic longitude
    deg : bool, optional
          degrees input/output  (False: radians in/out)

    Results
    -------

    (Vector)

    n : "ndarray"
        North NED coordinate (meters)
    e : "ndarray"
        East NED coordinate (meters)
    d : "ndarray"
        Down NED coordinate (meters)
    """
    e, n, u = ecef2enuv(x, y, z, lat0, lon0, deg=deg)

    return n, e, -u
