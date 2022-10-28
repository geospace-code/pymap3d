""" Transforms involving NED North East Down """
from __future__ import annotations

from typing import Any, Sequence, overload

try:
    from numpy import asarray
    from numpy.typing import NDArray
except ImportError:
    pass

from ._types import ArrayLike
from .ecef import ecef2enu, ecef2enuv, ecef2geodetic, enu2ecef
from .ellipsoid import Ellipsoid
from .enu import aer2enu, enu2aer, geodetic2enu


@overload
def aer2ned(
    az: float,
    elev: float,
    slantRange: float,
    deg: bool = True,
) -> tuple[float, float, float]:
    pass


@overload
def aer2ned(
    az: ArrayLike,
    elev: ArrayLike,
    slantRange: ArrayLike,
    deg: bool = True,
) -> tuple[NDArray[Any], NDArray[Any], NDArray[Any]]:
    pass


def aer2ned(
    az: float | ArrayLike, elev: float | ArrayLike, slantRange: float | ArrayLike, deg: bool = True
) -> tuple[float, float, float] | tuple[NDArray[Any], NDArray[Any], NDArray[Any]]:
    """
    converts azimuth, elevation, range to target from observer to North, East, Down

    Parameters
    -----------

    az : float
         azimuth
    elev : float
         elevation
    slantRange : float
         slant range [meters]
    deg : bool, optional
          degrees input/output  (False: radians in/out)

    Results
    -------
    n : float
        North NED coordinate (meters)
    e : float
        East NED coordinate (meters)
    d : float
        Down NED coordinate (meters)
    """
    e, n, u = aer2enu(az, elev, slantRange, deg=deg)  # type: ignore[arg-type]

    return n, e, -u


@overload
def ned2aer(
    n: float,
    e: float,
    d: float,
    deg: bool = True,
) -> tuple[float, float, float]:
    pass


@overload
def ned2aer(
    n: ArrayLike,
    e: ArrayLike,
    d: ArrayLike,
    deg: bool = True,
) -> tuple[NDArray[Any], NDArray[Any], NDArray[Any]]:
    pass


def ned2aer(
    n: float | ArrayLike, e: float | ArrayLike, d: float | ArrayLike, deg: bool = True
) -> tuple[float, float, float] | tuple[NDArray[Any], NDArray[Any], NDArray[Any]]:
    """
    converts North, East, Down to azimuth, elevation, range

    Parameters
    ----------

    n : float
        North NED coordinate (meters)
    e : float
        East NED coordinate (meters)
    d : float
        Down NED coordinate (meters)
    deg : bool, optional
        degrees input/output  (False: radians in/out)

    Results
    -------

    az : float
         azimuth
    elev : float
         elevation
    slantRange : float
         slant range [meters]
    """
    try:
        d = asarray(d)
    except NameError:
        pass
    assert not isinstance(d, Sequence)

    return enu2aer(e, n, -d, deg=deg)  # type: ignore[arg-type]


@overload
def ned2geodetic(
    n: float,
    e: float,
    d: float,
    lat0: float,
    lon0: float,
    h0: float,
    ell: Ellipsoid | None = None,
    deg: bool = True,
) -> tuple[float, float, float]:
    pass


@overload
def ned2geodetic(
    n: ArrayLike,
    e: ArrayLike,
    d: ArrayLike,
    lat0: ArrayLike,
    lon0: ArrayLike,
    h0: ArrayLike,
    ell: Ellipsoid | None = None,
    deg: bool = True,
) -> tuple[NDArray[Any], NDArray[Any], NDArray[Any]]:
    pass


def ned2geodetic(
    n: float | ArrayLike,
    e: float | ArrayLike,
    d: float | ArrayLike,
    lat0: float | ArrayLike,
    lon0: float | ArrayLike,
    h0: float | ArrayLike,
    ell: Ellipsoid | None = None,
    deg: bool = True,
) -> tuple[float, float, float] | tuple[NDArray[Any], NDArray[Any], NDArray[Any]]:
    """
    Converts North, East, Down to target latitude, longitude, altitude

    Parameters
    ----------

    n : float
        North NED coordinate (meters)
    e : float
        East NED coordinate (meters)
    d : float
        Down NED coordinate (meters)
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
        target geodetic latitude
    lon : float
        target geodetic longitude
    h : float
        target altitude above geodetic ellipsoid (meters)

    """
    try:
        d = asarray(d)
    except NameError:
        pass
    assert not isinstance(d, Sequence)

    x, y, z = enu2ecef(e, n, -d, lat0, lon0, h0, ell, deg=deg)  # type: ignore[misc, arg-type]

    return ecef2geodetic(x, y, z, ell, deg=deg)


@overload
def ned2ecef(
    n: float,
    e: float,
    d: float,
    lat0: float,
    lon0: float,
    h0: float,
    ell: Ellipsoid | None = None,
    deg: bool = True,
) -> tuple[float, float, float]:
    pass


@overload
def ned2ecef(
    n: ArrayLike,
    e: ArrayLike,
    d: ArrayLike,
    lat0: ArrayLike,
    lon0: ArrayLike,
    h0: ArrayLike,
    ell: Ellipsoid | None = None,
    deg: bool = True,
) -> tuple[NDArray[Any], NDArray[Any], NDArray[Any]]:
    pass


def ned2ecef(
    n: float | ArrayLike,
    e: float | ArrayLike,
    d: float | ArrayLike,
    lat0: float | ArrayLike,
    lon0: float | ArrayLike,
    h0: float | ArrayLike,
    ell: Ellipsoid | None = None,
    deg: bool = True,
) -> tuple[float, float, float] | tuple[NDArray[Any], NDArray[Any], NDArray[Any]]:
    """
    North, East, Down to target ECEF coordinates

    Parameters
    ----------

    n : float
        North NED coordinate (meters)
    e : float
        East NED coordinate (meters)
    d : float
        Down NED coordinate (meters)
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

    x : float
        ECEF x coordinate (meters)
    y : float
        ECEF y coordinate (meters)
    z : float
        ECEF z coordinate (meters)
    """
    try:
        d = asarray(d)
    except NameError:
        pass
    assert not isinstance(d, Sequence)

    return enu2ecef(e, n, -d, lat0, lon0, h0, ell, deg=deg)  # type: ignore[misc, arg-type]


@overload
def ecef2ned(
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
def ecef2ned(
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


def ecef2ned(
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
    Convert ECEF x,y,z to North, East, Down

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

    Results
    -------

    n : float
        North NED coordinate (meters)
    e : float
        East NED coordinate (meters)
    d : float
        Down NED coordinate (meters)

    """
    e, n, u = ecef2enu(x, y, z, lat0, lon0, h0, ell, deg=deg)  # type: ignore[misc, arg-type]

    return n, e, -u


@overload
def geodetic2ned(
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
def geodetic2ned(
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


def geodetic2ned(
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
    convert latitude, longitude, altitude of target to North, East, Down from observer

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


    Results
    -------

    n : float
        North NED coordinate (meters)
    e : float
        East NED coordinate (meters)
    d : float
        Down NED coordinate (meters)
    """
    e, n, u = geodetic2enu(lat, lon, h, lat0, lon0, h0, ell, deg=deg)  # type: ignore[misc, arg-type]

    return n, e, -u


@overload
def ecef2nedv(
    x: float,
    y: float,
    z: float,
    lat0: float,
    lon0: float,
    deg: bool = True,
) -> tuple[float, float, float]:
    pass


@overload
def ecef2nedv(
    x: ArrayLike,
    y: ArrayLike,
    z: ArrayLike,
    lat0: ArrayLike,
    lon0: ArrayLike,
    deg: bool = True,
) -> tuple[NDArray[Any], NDArray[Any], NDArray[Any]]:
    pass


def ecef2nedv(
    x: float | ArrayLike,
    y: float | ArrayLike,
    z: float | ArrayLike,
    lat0: float | ArrayLike,
    lon0: float | ArrayLike,
    deg: bool = True,
) -> tuple[float, float, float] | tuple[NDArray[Any], NDArray[Any], NDArray[Any]]:
    """
    for VECTOR between two points

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
    deg : bool, optional
          degrees input/output  (False: radians in/out)

    Results
    -------

    (Vector)

    n : float
        North NED coordinate (meters)
    e : float
        East NED coordinate (meters)
    d : float
        Down NED coordinate (meters)
    """
    e, n, u = ecef2enuv(x, y, z, lat0, lon0, deg=deg)  # type: ignore[misc, arg-type]

    return n, e, -u
