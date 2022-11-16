""" transforms involving AER: azimuth, elevation, slant range"""

from __future__ import annotations

from datetime import datetime
from typing import Any, overload

try:
    from numpy.typing import NDArray

    from .eci import ecef2eci, eci2ecef
except ImportError:
    pass

from ._types import ArrayLike
from .ecef import ecef2enu, ecef2geodetic, enu2uvw, geodetic2ecef
from .ellipsoid import Ellipsoid
from .enu import aer2enu, enu2aer, geodetic2enu

__all__ = ["aer2ecef", "ecef2aer", "geodetic2aer", "aer2geodetic", "eci2aer", "aer2eci"]


@overload
def ecef2aer(
    x: float,
    y: float,
    z: float,
    lat0: float,
    lon0: float,
    h0: float,
    ell: Ellipsoid | None = None,
    deg: bool = True,
) -> tuple[float, float, float]:
    pass


@overload
def ecef2aer(
    x: ArrayLike,
    y: ArrayLike,
    z: ArrayLike,
    lat0: ArrayLike,
    lon0: ArrayLike,
    h0: ArrayLike,
    ell: Ellipsoid | None = None,
    deg: bool = True,
) -> tuple[NDArray[Any], NDArray[Any], NDArray[Any]]:
    pass


def ecef2aer(
    x: float | ArrayLike,
    y: float | ArrayLike,
    z: float | ArrayLike,
    lat0: float | ArrayLike,
    lon0: float | ArrayLike,
    h0: float | ArrayLike,
    ell: Ellipsoid | None = None,
    deg: bool = True,
) -> tuple[float, float, float] | tuple[NDArray[Any], NDArray[Any], NDArray[Any]]:
    """
    compute azimuth, elevation and slant range from an Observer to a Point with ECEF coordinates.

    ECEF input location is with units of meters

    Parameters
    ----------

    x : float
        ECEF x coordinate (meters)
    y : float
        ECEF y coordinate (meters)
    z : float
        ECEF z coordinate (meters)
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

    Returns
    -------
    az : float
         azimuth to target
    el : float
         elevation to target
    srange : float
         slant range [meters]
    """
    xEast, yNorth, zUp = ecef2enu(x, y, z, lat0, lon0, h0, ell, deg=deg)  # type: ignore[misc, arg-type]

    return enu2aer(xEast, yNorth, zUp, deg=deg)


@overload
def geodetic2aer(
    lat: float,
    lon: float,
    h: float,
    lat0: float,
    lon0: float,
    h0: float,
    ell: Ellipsoid | None = None,
    deg: bool = True,
) -> tuple[float, float, float]:
    pass


@overload
def geodetic2aer(
    lat: ArrayLike,
    lon: ArrayLike,
    h: ArrayLike,
    lat0: ArrayLike,
    lon0: ArrayLike,
    h0: ArrayLike,
    ell: Ellipsoid | None = None,
    deg: bool = True,
) -> tuple[NDArray[Any], NDArray[Any], NDArray[Any]]:
    pass


def geodetic2aer(
    lat: float | ArrayLike,
    lon: float | ArrayLike,
    h: float | ArrayLike,
    lat0: float | ArrayLike,
    lon0: float | ArrayLike,
    h0: float | ArrayLike,
    ell: Ellipsoid | None = None,
    deg: bool = True,
) -> tuple[float, float, float] | tuple[NDArray[Any], NDArray[Any], NDArray[Any]]:
    """
    gives azimuth, elevation and slant range from an Observer to a Point with geodetic coordinates.


    Parameters
    ----------

    lat : float
        target geodetic latitude
    lon : float
        target geodetic longitude
    h : float
        target altitude above geodetic ellipsoid (meters)
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

    Returns
    -------
    az : float
         azimuth
    el : float
         elevation
    srange : float
         slant range [meters]
    """
    e, n, u = geodetic2enu(lat, lon, h, lat0, lon0, h0, ell, deg=deg)  # type: ignore[misc, arg-type]

    return enu2aer(e, n, u, deg=deg)


@overload
def aer2geodetic(
    az: float,
    el: float,
    srange: float,
    lat0: float,
    lon0: float,
    h0: float,
    ell: Ellipsoid | None = None,
    deg: bool = True,
) -> tuple[float, float, float]:
    pass


@overload
def aer2geodetic(
    az: ArrayLike,
    el: ArrayLike,
    srange: ArrayLike,
    lat0: ArrayLike,
    lon0: ArrayLike,
    h0: ArrayLike,
    ell: Ellipsoid | None = None,
    deg: bool = True,
) -> tuple[NDArray[Any], NDArray[Any], NDArray[Any]]:
    pass


def aer2geodetic(
    az: float | ArrayLike,
    el: float | ArrayLike,
    srange: float | ArrayLike,
    lat0: float | ArrayLike,
    lon0: float | ArrayLike,
    h0: float | ArrayLike,
    ell: Ellipsoid | None = None,
    deg: bool = True,
) -> tuple[float, float, float] | tuple[NDArray[Any], NDArray[Any], NDArray[Any]]:
    """
    gives geodetic coordinates of a point with az, el, range
    from an observer at lat0, lon0, h0

    Parameters
    ----------
    az : float
         azimuth to target
    el : float
         elevation to target
    srange : float
         slant range [meters]
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

    Returns
    -------

    In reference ellipsoid system:

    lat : float
          geodetic latitude
    lon : float
          geodetic longitude
    alt : float
          altitude above ellipsoid  (meters)
    """
    x, y, z = aer2ecef(az, el, srange, lat0, lon0, h0, ell=ell, deg=deg)  # type: ignore[misc, arg-type]

    return ecef2geodetic(x, y, z, ell=ell, deg=deg)


@overload
def eci2aer(
    x: float,
    y: float,
    z: float,
    lat0: float,
    lon0: float,
    h0: float,
    t: datetime,
    *,
    deg: bool = True,
) -> tuple[float, float, float]:
    pass


@overload
def eci2aer(
    x: ArrayLike,
    y: ArrayLike,
    z: ArrayLike,
    lat0: ArrayLike,
    lon0: ArrayLike,
    h0: ArrayLike,
    t: datetime,
    *,
    deg: bool = True,
) -> tuple[NDArray[Any], NDArray[Any], NDArray[Any]]:
    pass


def eci2aer(
    x: float | ArrayLike,
    y: float | ArrayLike,
    z: float | ArrayLike,
    lat0: float | ArrayLike,
    lon0: float | ArrayLike,
    h0: float | ArrayLike,
    t: datetime,
    *,
    deg: bool = True,
) -> tuple[float, float, float] | tuple[NDArray[Any], NDArray[Any], NDArray[Any]]:
    """
    takes Earth Centered Inertial x,y,z ECI coordinates of point and gives az, el, slant range from Observer

    Parameters
    ----------

    x : float
        ECI x-location [meters]
    y : float
        ECI y-location [meters]
    z : float
        ECI z-location [meters]
    lat0 : float
           Observer geodetic latitude
    lon0 : float
           Observer geodetic longitude
    h0 : float
         observer altitude above geodetic ellipsoid (meters)
    t : datetime.datetime
        Observation time
    deg : bool, optional
        true: degrees, false: radians

    Returns
    -------
    az : float
         azimuth to target
    el : float
         elevation to target
    srange : float
         slant range [meters]
    """

    try:
        xecef, yecef, zecef = eci2ecef(x, y, z, t)  # type: ignore[arg-type]
    except NameError:
        raise ImportError("pip install numpy")

    return ecef2aer(xecef, yecef, zecef, lat0, lon0, h0, deg=deg)  # type: ignore[arg-type]


@overload
def aer2eci(
    az: float,
    el: float,
    srange: float,
    lat0: float,
    lon0: float,
    h0: float,
    t: datetime,
    ell: Ellipsoid | None = None,
    *,
    deg: bool = True,
) -> tuple[float, float, float]:
    pass


@overload
def aer2eci(
    az: ArrayLike,
    el: ArrayLike,
    srange: ArrayLike,
    lat0: ArrayLike,
    lon0: ArrayLike,
    h0: ArrayLike,
    t: datetime,
    ell: Ellipsoid | None = None,
    *,
    deg: bool = True,
) -> tuple[NDArray[Any], NDArray[Any], NDArray[Any]]:
    pass


def aer2eci(
    az: float | ArrayLike,
    el: float | ArrayLike,
    srange: float | ArrayLike,
    lat0: float | ArrayLike,
    lon0: float | ArrayLike,
    h0: float | ArrayLike,
    t: datetime,
    ell: Ellipsoid | None = None,
    *,
    deg: bool = True,
) -> tuple[float, float, float] | tuple[NDArray[Any], NDArray[Any], NDArray[Any]]:
    """
    gives ECI of a point from an observer at az, el, slant range

    Parameters
    ----------
    az : float
         azimuth to target
    el : float
         elevation to target
    srange : float
         slant range [meters]
    lat0 : float
           Observer geodetic latitude
    lon0 : float
           Observer geodetic longitude
    h0 : float
         observer altitude above geodetic ellipsoid (meters)
    t : datetime.datetime
        Observation time
    ell : Ellipsoid, optional
          reference ellipsoid
    deg : bool, optional
          degrees input/output  (False: radians in/out)

    Returns
    -------

    Earth Centered Inertial x,y,z

    x : float
        ECEF x coordinate (meters)
    y : float
        ECEF y coordinate (meters)
    z : float
        ECEF z coordinate (meters)
    """

    x, y, z = aer2ecef(az, el, srange, lat0, lon0, h0, ell, deg=deg)  # type: ignore[misc, arg-type]

    try:
        return ecef2eci(x, y, z, t)
    except NameError:
        raise ImportError("pip install numpy")


@overload
def aer2ecef(
    az: float,
    el: float,
    srange: float,
    lat0: float,
    lon0: float,
    alt0: float,
    ell: Ellipsoid | None = None,
    deg: bool = True,
) -> tuple[float, float, float]:
    pass


@overload
def aer2ecef(
    az: ArrayLike,
    el: ArrayLike,
    srange: ArrayLike,
    lat0: ArrayLike,
    lon0: ArrayLike,
    alt0: ArrayLike,
    ell: Ellipsoid | None = None,
    deg: bool = True,
) -> tuple[NDArray[Any], NDArray[Any], NDArray[Any]]:
    pass


def aer2ecef(
    az: float | ArrayLike,
    el: float | ArrayLike,
    srange: float | ArrayLike,
    lat0: float | ArrayLike,
    lon0: float | ArrayLike,
    alt0: float | ArrayLike,
    ell: Ellipsoid | None = None,
    deg: bool = True,
) -> tuple[float, float, float] | tuple[NDArray[Any], NDArray[Any], NDArray[Any]]:
    """
    converts target azimuth, elevation, range from observer at lat0,lon0,alt0 to ECEF coordinates.

    Parameters
    ----------
    az : float
         azimuth to target
    el : float
         elevation to target
    srange : float
         slant range [meters]
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

    Returns
    -------

    ECEF (Earth centered, Earth fixed)  x,y,z

    x : float
        ECEF x coordinate (meters)
    y : float
        ECEF y coordinate (meters)
    z : float
        ECEF z coordinate (meters)


    Notes
    ------
    if srange==NaN, z=NaN
    """
    # Origin of the local system in geocentric coordinates.
    x0, y0, z0 = geodetic2ecef(lat0, lon0, alt0, ell, deg=deg)  # type: ignore[misc, arg-type]
    # Convert Local Spherical AER to ENU
    e1, n1, u1 = aer2enu(az, el, srange, deg=deg)  # type: ignore[arg-type]
    # Rotating ENU to ECEF
    dx, dy, dz = enu2uvw(e1, n1, u1, lat0, lon0, deg=deg)  # type: ignore[arg-type]
    # Origin + offset from origin equals position in ECEF
    return x0 + dx, y0 + dy, z0 + dz
