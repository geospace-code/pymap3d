"""transforms involving AER: azimuth, elevation, slant range"""

from __future__ import annotations

from datetime import datetime

from .ecef import ecef2enu, ecef2geodetic, enu2uvw, geodetic2ecef
from .ellipsoid import Ellipsoid
from .enu import aer2enu, enu2aer, geodetic2enu
from ._typing import FloatLike

try:
    from .eci import ecef2eci, eci2ecef
except ImportError:
    from ._dummy import ecef2eci, eci2ecef


__all__ = ["aer2ecef", "ecef2aer", "geodetic2aer", "aer2geodetic", "eci2aer", "aer2eci"]


def ecef2aer(
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
    compute azimuth, elevation and slant range from an Observer to a Point with ECEF coordinates.

    ECEF input location is with units of meters

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

    Returns
    -------
    az : array-like float
         azimuth to target
    el : array-like float
         elevation to target
    srange : array-like float
         slant range [meters]
    """

    xEast, yNorth, zUp = ecef2enu(x, y, z, lat0, lon0, h0, ell, deg=deg)

    return enu2aer(xEast, yNorth, zUp, deg=deg)


def geodetic2aer(
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
    gives azimuth, elevation and slant range from an Observer to a Point with geodetic coordinates.


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
         observer altitude above geodetic ellipsoid (meters)
    ell : Ellipsoid, optional
          reference ellipsoid
    deg : bool, optional
          degrees input/output  (False: radians in/out)

    Returns
    -------
    az : array-like float
         azimuth
    el : array-like float
         elevation
    srange : array-like float
         slant range [meters]
    """
    e, n, u = geodetic2enu(lat, lon, h, lat0, lon0, h0, ell, deg=deg)

    return enu2aer(e, n, u, deg=deg)


def aer2geodetic(
    az: FloatLike,
    el: FloatLike,
    srange: FloatLike,
    lat0: FloatLike,
    lon0: FloatLike,
    h0: FloatLike,
    ell: Ellipsoid | None = None,
    deg: bool = True,
) -> tuple:
    """
    gives geodetic coordinates of a point with az, el, range
    from an observer at lat0, lon0, h0

    Parameters
    ----------
    az : array-like float
         azimuth to target
    el : array-like float
         elevation to target
    srange : array-like float
         slant range [meters]
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

    Returns
    -------

    In reference ellipsoid system:

    lat : array-like float
          geodetic latitude
    lon : array-like float
          geodetic longitude
    alt : array-like float
          altitude above ellipsoid  (meters)
    """
    x, y, z = aer2ecef(az, el, srange, lat0, lon0, h0, ell=ell, deg=deg)

    return ecef2geodetic(x, y, z, ell=ell, deg=deg)


def eci2aer(
    x: FloatLike,
    y: FloatLike,
    z: FloatLike,
    lat0: FloatLike,
    lon0: FloatLike,
    h0: FloatLike,
    t: datetime,
    *,
    deg: bool = True,
    delta_ut1: float = 0.0,
    xp: float = 0.0,
    yp: float = 0.0,
) -> tuple:
    """
    takes Earth Centered Inertial x,y,z ECI coordinates of point and gives az, el,
    slant range from Observer

    Parameters
    ----------

    x : array-like float
        ECI x-location [meters]
    y : array-like float
        ECI y-location [meters]
    z : array-like float
        ECI z-location [meters]
    lat0 : array-like float
           Observer geodetic latitude
    lon0 : array-like float
           Observer geodetic longitude
    h0 : array-like float
         observer altitude above geodetic ellipsoid (meters)
    t : datetime.datetime
        Observation time
    deg : bool, optional
        true: degrees, false: radians
    delta_ut1 : float, optional
        UT1-UTC in seconds for the pure-Python path. Defaults to ``0.0``.
    xp : float, optional
        Polar motion x coordinate in arcseconds for the pure-Python path.
    yp : float, optional
        Polar motion y coordinate in arcseconds for the pure-Python path.

    Returns
    -------
    az : array-like float
         azimuth to target
    el : array-like float
         elevation to target
    srange : array-like float
         slant range [meters]
    """

    xe, ye, ze = eci2ecef(x, y, z, t, delta_ut1=delta_ut1, xp=xp, yp=yp)

    return ecef2aer(xe, ye, ze, lat0, lon0, h0, deg=deg)


def aer2eci(
    az: FloatLike,
    el: FloatLike,
    srange: FloatLike,
    lat0: FloatLike,
    lon0: FloatLike,
    h0: FloatLike,
    t: datetime,
    ell: Ellipsoid | None = None,
    *,
    deg: bool = True,
    delta_ut1: float = 0.0,
    xp: float = 0.0,
    yp: float = 0.0,
) -> tuple:
    """
    gives ECI of a point from an observer at az, el, slant range

    Parameters
    ----------
    az : array-like float
         azimuth to target
    el : array-like float
         elevation to target
    srange : array-like float
         slant range [meters]
    lat0 : array-like float
           Observer geodetic latitude
    lon0 : array-like float
           Observer geodetic longitude
    h0 : array-like float
         observer altitude above geodetic ellipsoid (meters)
    t : datetime.datetime
        Observation time
    ell : Ellipsoid, optional
          reference ellipsoid
    deg : bool, optional
          degrees input/output  (False: radians in/out)
    delta_ut1 : float, optional
        UT1-UTC in seconds for the pure-Python path. Defaults to ``0.0``.
    xp : float, optional
        Polar motion x coordinate in arcseconds for the pure-Python path.
    yp : float, optional
        Polar motion y coordinate in arcseconds for the pure-Python path.

    Returns
    -------

    Earth Centered Inertial x,y,z

    x : array-like float
        ECEF x coordinate (meters)
    y : array-like float
        ECEF y coordinate (meters)
    z : array-like float
        ECEF z coordinate (meters)
    """

    x, y, z = aer2ecef(az, el, srange, lat0, lon0, h0, ell, deg=deg)

    return ecef2eci(x, y, z, t, delta_ut1=delta_ut1, xp=xp, yp=yp)


def aer2ecef(
    az: FloatLike,
    el: FloatLike,
    srange: FloatLike,
    lat0: FloatLike,
    lon0: FloatLike,
    alt0: FloatLike,
    ell: Ellipsoid | None = None,
    deg: bool = True,
) -> tuple:
    """
    converts target azimuth, elevation, range from observer at lat0,lon0,alt0 to ECEF coordinates.

    Parameters
    ----------
    az : array-like float
         azimuth to target
    el : array-like float
         elevation to target
    srange : array-like float
         slant range [meters]
    lat0 : array-like float
           Observer geodetic latitude
    lon0 : array-like float
           Observer geodetic longitude
    alt0 : array-like float
         observer altitude above geodetic ellipsoid (meters)
    ell : Ellipsoid, optional
          reference ellipsoid
    deg : bool, optional
          degrees input/output  (False: radians in/out)

    Returns
    -------

    ECEF (Earth centered, Earth fixed)  x,y,z

    x : array-like float
        ECEF x coordinate (meters)
    y : array-like float
        ECEF y coordinate (meters)
    z : array-like float
        ECEF z coordinate (meters)


    Notes
    ------
    if srange==NaN, z=NaN
    """
    # Origin of the local system in geocentric coordinates.
    x0, y0, z0 = geodetic2ecef(lat0, lon0, alt0, ell, deg=deg)
    # Convert Local Spherical AER to ENU
    e1, n1, u1 = aer2enu(az, el, srange, deg=deg)
    # Rotating ENU to ECEF
    dx, dy, dz = enu2uvw(e1, n1, u1, lat0, lon0, deg=deg)
    # Origin + offset from origin equals position in ECEF
    return x0 + dx, y0 + dy, z0 + dz
