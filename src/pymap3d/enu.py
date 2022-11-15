""" transforms involving ENU East North Up """
from __future__ import annotations

from math import tau
from typing import Any, overload

try:
    from numpy import asarray
    from numpy.typing import NDArray
except ImportError:
    pass
from ._types import ArrayLike
from .ecef import ecef2geodetic, enu2ecef, geodetic2ecef, uvw2enu
from .ellipsoid import Ellipsoid
from .mathfun import atan2, cos, degrees, hypot, radians, sin

__all__ = ["enu2aer", "aer2enu", "enu2geodetic", "geodetic2enu"]


@overload
def enu2aer(
    e: float,
    n: float,
    u: float,
    deg: bool = True,
) -> tuple[float, float, float]:
    pass


@overload
def enu2aer(
    e: ArrayLike,
    n: ArrayLike,
    u: ArrayLike,
    deg: bool = True,
) -> tuple[NDArray[Any], NDArray[Any], NDArray[Any]]:
    pass


def enu2aer(
    e: float | ArrayLike, n: float | ArrayLike, u: float | ArrayLike, deg: bool = True
) -> tuple[float, float, float] | tuple[NDArray[Any], NDArray[Any], NDArray[Any]]:
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
        e = asarray(e)
        n = asarray(n)
        u = asarray(u)
        e[abs(e) < 1e-3] = 0.0
        n[abs(n) < 1e-3] = 0.0
        u[abs(u) < 1e-3] = 0.0
    except (TypeError, NameError):
        assert isinstance(e, float) and isinstance(n, float) and isinstance(u, float)
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


@overload
def aer2enu(
    az: float,
    el: float,
    srange: float,
    deg: bool = True,
) -> tuple[float, float, float]:
    pass


@overload
def aer2enu(
    az: ArrayLike,
    el: ArrayLike,
    srange: ArrayLike,
    deg: bool = True,
) -> tuple[NDArray[Any], NDArray[Any], NDArray[Any]]:
    pass


def aer2enu(
    az: float | ArrayLike, el: float | ArrayLike, srange: float | ArrayLike, deg: bool = True
) -> tuple[float, float, float] | tuple[NDArray[Any], NDArray[Any], NDArray[Any]]:
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
        assert isinstance(srange, int) or isinstance(srange, float)
        if srange < 0:
            raise ValueError("Slant range  [0, Infinity)")

    r = srange * cos(el)

    return r * sin(az), r * cos(az), srange * sin(el)


@overload
def enu2geodetic(
    e: float,
    n: float,
    u: float,
    lat0: float,
    lon0: float,
    h0: float,
    ell: Ellipsoid | None = None,
    deg: bool = True,
) -> tuple[float, float, float]:
    pass


@overload
def enu2geodetic(
    e: ArrayLike,
    n: ArrayLike,
    u: ArrayLike,
    lat0: ArrayLike,
    lon0: ArrayLike,
    h0: ArrayLike,
    ell: Ellipsoid | None = None,
    deg: bool = True,
) -> tuple[NDArray[Any], NDArray[Any], NDArray[Any]]:
    pass


def enu2geodetic(
    e: float | ArrayLike,
    n: float | ArrayLike,
    u: float | ArrayLike,
    lat0: float | ArrayLike,
    lon0: float | ArrayLike,
    h0: float | ArrayLike,
    ell: Ellipsoid | None = None,
    deg: bool = True,
) -> tuple[float, float, float] | tuple[NDArray[Any], NDArray[Any], NDArray[Any]]:
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

    x, y, z = enu2ecef(e, n, u, lat0, lon0, h0, ell, deg=deg)  # type: ignore[misc, arg-type]

    return ecef2geodetic(x, y, z, ell, deg=deg)


@overload
def geodetic2enu(
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
def geodetic2enu(
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


def geodetic2enu(
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
    x1, y1, z1 = geodetic2ecef(lat, lon, h, ell, deg=deg)  # type: ignore[misc, arg-type]
    x2, y2, z2 = geodetic2ecef(lat0, lon0, h0, ell, deg=deg)  # type: ignore[misc, arg-type]

    return uvw2enu(x1 - x2, y1 - y2, z1 - z2, lat0, lon0, deg=deg)  # type: ignore[arg-type]
