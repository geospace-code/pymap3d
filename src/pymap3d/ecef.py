""" Transforms involving ECEF: earth-centered, earth-fixed frame """
from __future__ import annotations

from typing import Any, Sequence, overload

try:
    from numpy import asarray, finfo, where
    from numpy.typing import NDArray

    from .eci import ecef2eci, eci2ecef
except ImportError:
    pass

from datetime import datetime
from math import pi

from ._types import ArrayLike
from .ellipsoid import Ellipsoid
from .mathfun import atan, atan2, cos, degrees, hypot, radians, sin, sqrt, tan
from .utils import sanitize

__all__ = [
    "geodetic2ecef",
    "ecef2geodetic",
    "ecef2enuv",
    "ecef2enu",
    "enu2uvw",
    "uvw2enu",
    "eci2geodetic",
    "geodetic2eci",
    "enu2ecef",
]


@overload
def geodetic2ecef(
    lat: float,
    lon: float,
    alt: float,
    ell: Ellipsoid | None = None,
    deg: bool = True,
) -> tuple[float, float, float]:
    pass


@overload
def geodetic2ecef(
    lat: ArrayLike,
    lon: ArrayLike,
    alt: ArrayLike,
    ell: Ellipsoid | None = None,
    deg: bool = True,
) -> tuple[NDArray[Any], NDArray[Any], NDArray[Any]]:
    pass


def geodetic2ecef(
    lat: float | ArrayLike,
    lon: float | ArrayLike,
    alt: float | ArrayLike,
    ell: Ellipsoid | None = None,
    deg: bool = True,
) -> tuple[float, float, float] | tuple[NDArray[Any], NDArray[Any], NDArray[Any]]:
    """
    point transformation from Geodetic of specified ellipsoid (default WGS-84) to ECEF

    Parameters
    ----------

    lat
           target geodetic latitude
    lon
           target geodetic longitude
    h
         target altitude above geodetic ellipsoid (meters)
    ell : Ellipsoid, optional
          reference ellipsoid
    deg : bool, optional
          degrees input/output  (False: radians in/out)


    Returns
    -------

    ECEF (Earth centered, Earth fixed)  x,y,z

    x
        target x ECEF coordinate (meters)
    y
        target y ECEF coordinate (meters)
    z
        target z ECEF coordinate (meters)
    """
    lat, ell = sanitize(lat, ell, deg)
    if deg:
        lon = radians(lon)

    # radius of curvature of the prime vertical section
    N = ell.semimajor_axis**2 / sqrt(
        ell.semimajor_axis**2 * cos(lat) ** 2 + ell.semiminor_axis**2 * sin(lat) ** 2
    )
    # Compute cartesian (geocentric) coordinates given (curvilinear) geodetic coordinates.
    x = (N + alt) * cos(lat) * cos(lon)
    y = (N + alt) * cos(lat) * sin(lon)
    z = (N * (ell.semiminor_axis / ell.semimajor_axis) ** 2 + alt) * sin(lat)

    return x, y, z


@overload
def ecef2geodetic(
    x: float,
    y: float,
    z: float,
    ell: Ellipsoid | None = None,
    deg: bool = True,
) -> tuple[float, float, float]:
    pass


@overload
def ecef2geodetic(
    x: ArrayLike,
    y: ArrayLike,
    z: ArrayLike,
    ell: Ellipsoid | None = None,
    deg: bool = True,
) -> tuple[NDArray[Any], NDArray[Any], NDArray[Any]]:
    pass


def ecef2geodetic(
    x: float | ArrayLike,
    y: float | ArrayLike,
    z: float | ArrayLike,
    ell: Ellipsoid | None = None,
    deg: bool = True,
) -> tuple[float, float, float] | tuple[NDArray[Any], NDArray[Any], NDArray[Any]]:
    """
    convert ECEF (meters) to geodetic coordinates

    Parameters
    ----------
    x
        target x ECEF coordinate (meters)
    y
        target y ECEF coordinate (meters)
    z
        target z ECEF coordinate (meters)
    ell : Ellipsoid, optional
          reference ellipsoid
    deg : bool, optional
          degrees input/output  (False: radians in/out)

    Returns
    -------
    lat
           target geodetic latitude
    lon
           target geodetic longitude
    alt
         target altitude above geodetic ellipsoid (meters)

    based on:
    You, Rey-Jer. (2000). Transformation of Cartesian to Geodetic Coordinates without Iterations.
    Journal of Surveying Engineering. doi: 10.1061/(ASCE)0733-9453
    """

    if ell is None:
        ell = Ellipsoid.from_name("wgs84")

    try:
        x = asarray(x)
        y = asarray(y)
        z = asarray(z)
    except NameError:
        pass

    assert (
        not isinstance(x, Sequence) and not isinstance(y, Sequence) and not isinstance(z, Sequence)
    )

    r = sqrt(x**2 + y**2 + z**2)

    E = sqrt(ell.semimajor_axis**2 - ell.semiminor_axis**2)

    # eqn. 4a
    u = sqrt(0.5 * (r**2 - E**2) + 0.5 * sqrt((r**2 - E**2) ** 2 + 4 * E**2 * z**2))

    Q = hypot(x, y)

    huE = hypot(u, E)

    # eqn. 4b
    try:
        Beta = atan(huE / u * z / hypot(x, y))
    except ZeroDivisionError:
        if z >= 0:
            Beta = pi / 2
        else:
            Beta = -pi / 2

    # eqn. 13
    dBeta = ((ell.semiminor_axis * u - ell.semimajor_axis * huE + E**2) * sin(Beta)) / (
        ell.semimajor_axis * huE * 1 / cos(Beta) - E**2 * cos(Beta)
    )

    Beta += dBeta

    # eqn. 4c
    # %% final output
    lat = atan(ell.semimajor_axis / ell.semiminor_axis * tan(Beta))

    try:
        # patch latitude for float32 precision loss
        lim_pi2 = pi / 2 - finfo(dBeta.dtype).eps
        lat = where(Beta >= lim_pi2, pi / 2, lat)
        lat = where(Beta <= -lim_pi2, -pi / 2, lat)
    except NameError:
        pass

    lon = atan2(y, x)

    # eqn. 7
    cosBeta = cos(Beta)
    try:
        # patch altitude for float32 precision loss
        cosBeta = where(Beta >= lim_pi2, 0, cosBeta)
        cosBeta = where(Beta <= -lim_pi2, 0, cosBeta)
    except NameError:
        pass

    alt = hypot(z - ell.semiminor_axis * sin(Beta), Q - ell.semimajor_axis * cosBeta)

    # inside ellipsoid?
    inside = (
        x**2 / ell.semimajor_axis**2
        + y**2 / ell.semimajor_axis**2
        + z**2 / ell.semiminor_axis**2
        < 1
    )

    try:
        if inside.any():  # type: ignore[union-attr]
            # avoid all false assignment bug
            alt[inside] = -alt[inside]
    except (TypeError, AttributeError):
        if inside:
            alt = -alt

    if deg:
        lat = degrees(lat)
        lon = degrees(lon)

    return lat, lon, alt


@overload
def ecef2enuv(
    u: float,
    v: float,
    w: float,
    lat0: float,
    lon0: float,
    deg: bool = True,
) -> tuple[float, float, float]:
    pass


@overload
def ecef2enuv(
    u: ArrayLike,
    v: ArrayLike,
    w: ArrayLike,
    lat0: ArrayLike,
    lon0: ArrayLike,
    deg: bool = True,
) -> tuple[NDArray[Any], NDArray[Any], NDArray[Any]]:
    pass


def ecef2enuv(
    u: float | ArrayLike,
    v: float | ArrayLike,
    w: float | ArrayLike,
    lat0: float | ArrayLike,
    lon0: float | ArrayLike,
    deg: bool = True,
) -> tuple[float, float, float] | tuple[NDArray[Any], NDArray[Any], NDArray[Any]]:
    """
    VECTOR from observer to target  ECEF => ENU

    Parameters
    ----------
    u
        target x ECEF coordinate (meters)
    v
        target y ECEF coordinate (meters)
    w
        target z ECEF coordinate (meters)
    lat0
           Observer geodetic latitude
    lon0
           Observer geodetic longitude
    h0
         observer altitude above geodetic ellipsoid (meters)
    deg : bool, optional
          degrees input/output  (False: radians in/out)

    Returns
    -------
    uEast
        target east ENU coordinate (meters)
    vNorth
        target north ENU coordinate (meters)
    wUp
        target up ENU coordinate (meters)

    """
    if deg:
        lat0 = radians(lat0)
        lon0 = radians(lon0)

    t = cos(lon0) * u + sin(lon0) * v
    uEast = -sin(lon0) * u + cos(lon0) * v
    wUp = cos(lat0) * t + sin(lat0) * w
    vNorth = -sin(lat0) * t + cos(lat0) * w

    return uEast, vNorth, wUp


@overload
def ecef2enu(
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
def ecef2enu(
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


def ecef2enu(
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
    from observer to target, ECEF => ENU

    Parameters
    ----------
    x
        target x ECEF coordinate (meters)
    y
        target y ECEF coordinate (meters)
    z
        target z ECEF coordinate (meters)
    lat0
           Observer geodetic latitude
    lon0
           Observer geodetic longitude
    h0
         observer altitude above geodetic ellipsoid (meters)
    ell : Ellipsoid, optional
          reference ellipsoid
    deg : bool, optional
          degrees input/output  (False: radians in/out)

    Returns
    -------
    East
        target east ENU coordinate (meters)
    North
        target north ENU coordinate (meters)
    Up
        target up ENU coordinate (meters)

    """

    x0, y0, z0 = geodetic2ecef(lat0, lon0, h0, ell, deg=deg)  # type: ignore[misc, arg-type]

    return uvw2enu(x - x0, y - y0, z - z0, lat0, lon0, deg=deg)  # type: ignore[misc, arg-type, operator]


@overload
def enu2uvw(
    east: float,
    north: float,
    up: float,
    lat0: float,
    lon0: float,
    deg: bool = True,
) -> tuple[float, float, float]:
    pass


@overload
def enu2uvw(
    east: ArrayLike,
    north: ArrayLike,
    up: ArrayLike,
    lat0: ArrayLike,
    lon0: ArrayLike,
    deg: bool = True,
) -> tuple[NDArray[Any], NDArray[Any], NDArray[Any]]:
    pass


def enu2uvw(
    east: float | ArrayLike,
    north: float | ArrayLike,
    up: float | ArrayLike,
    lat0: float | ArrayLike,
    lon0: float | ArrayLike,
    deg: bool = True,
) -> tuple[float, float, float] | tuple[NDArray[Any], NDArray[Any], NDArray[Any]]:
    """
    Parameters
    ----------

    e1
        target east ENU coordinate (meters)
    n1
        target north ENU coordinate (meters)
    u1
        target up ENU coordinate (meters)

    Results
    -------

    u
    v
    w
    """

    if deg:
        lat0 = radians(lat0)
        lon0 = radians(lon0)

    t = cos(lat0) * up - sin(lat0) * north
    w = sin(lat0) * up + cos(lat0) * north

    u = cos(lon0) * t - sin(lon0) * east
    v = sin(lon0) * t + cos(lon0) * east

    return u, v, w


@overload
def uvw2enu(
    u: float,
    v: float,
    w: float,
    lat0: float,
    lon0: float,
    deg: bool = True,
) -> tuple[float, float, float]:
    pass


@overload
def uvw2enu(
    u: ArrayLike,
    v: ArrayLike,
    w: ArrayLike,
    lat0: ArrayLike,
    lon0: ArrayLike,
    deg: bool = True,
) -> tuple[NDArray[Any], NDArray[Any], NDArray[Any]]:
    pass


def uvw2enu(
    u: float | ArrayLike,
    v: float | ArrayLike,
    w: float | ArrayLike,
    lat0: float | ArrayLike,
    lon0: float | ArrayLike,
    deg: bool = True,
) -> tuple[float, float, float] | tuple[NDArray[Any], NDArray[Any], NDArray[Any]]:
    """
    Parameters
    ----------

    u
    v
    w


    Results
    -------

    East
        target east ENU coordinate (meters)
    North
        target north ENU coordinate (meters)
    Up
        target up ENU coordinate (meters)
    """
    if deg:
        lat0 = radians(lat0)
        lon0 = radians(lon0)

    t = cos(lon0) * u + sin(lon0) * v
    East = -sin(lon0) * u + cos(lon0) * v
    Up = cos(lat0) * t + sin(lat0) * w
    North = -sin(lat0) * t + cos(lat0) * w

    return East, North, Up


@overload
def eci2geodetic(
    x: float, y: float, z: float, t: datetime, ell: Ellipsoid | None = None, *, deg: bool = True
) -> tuple[float, float, float]:
    pass


@overload
def eci2geodetic(
    x: ArrayLike,
    y: ArrayLike,
    z: ArrayLike,
    t: datetime,
    ell: Ellipsoid | None = None,
    *,
    deg: bool = True,
) -> tuple[NDArray[Any], NDArray[Any], NDArray[Any]]:
    pass


def eci2geodetic(
    x: float | ArrayLike,
    y: float | ArrayLike,
    z: float | ArrayLike,
    t: datetime,
    ell: Ellipsoid | None = None,
    *,
    deg: bool = True,
) -> tuple[float, float, float] | tuple[NDArray[Any], NDArray[Any], NDArray[Any]]:
    """
    convert Earth Centered Internal ECI to geodetic coordinates

    J2000 time

    Parameters
    ----------
    x
        ECI x-location [meters]
    y
        ECI y-location [meters]
    z
        ECI z-location [meters]
    t : datetime.datetime, float
        UTC time
    ell : Ellipsoid, optional
        planet ellipsoid model
    deg : bool, optional
        if True, degrees. if False, radians

    Results
    -------
    lat
          geodetic latitude
    lon
          geodetic longitude
    alt
          altitude above ellipsoid  (meters)

    eci2geodetic() a.k.a. eci2lla()
    """

    try:
        xecef, yecef, zecef = eci2ecef(x, y, z, t)  # type: ignore[arg-type]
    except NameError:
        raise ImportError("pip install numpy")

    return ecef2geodetic(xecef, yecef, zecef, ell, deg)


@overload
def geodetic2eci(
    lat: float,
    lon: float,
    alt: float,
    t: datetime,
    ell: Ellipsoid | None = None,
    *,
    deg: bool = True,
) -> tuple[float, float, float]:
    pass


@overload
def geodetic2eci(
    lat: ArrayLike,
    lon: ArrayLike,
    alt: ArrayLike,
    t: datetime,
    ell: Ellipsoid | None = None,
    *,
    deg: bool = True,
) -> tuple[NDArray[Any], NDArray[Any], NDArray[Any]]:
    pass


def geodetic2eci(
    lat: float | ArrayLike,
    lon: float | ArrayLike,
    alt: float | ArrayLike,
    t: datetime,
    ell: Ellipsoid | None = None,
    *,
    deg: bool = True,
) -> tuple[float, float, float] | tuple[NDArray[Any], NDArray[Any], NDArray[Any]]:
    """
    convert geodetic coordinates to Earth Centered Internal ECI

    J2000 frame

    Parameters
    ----------
    lat
        geodetic latitude
    lon
        geodetic longitude
    alt
        altitude above ellipsoid  (meters)
    t : datetime.datetime, float
        UTC time
    ell : Ellipsoid, optional
        planet ellipsoid model
    deg : bool, optional
        if True, degrees. if False, radians

    Results
    -------
    x
        ECI x-location [meters]
    y
        ECI y-location [meters]
    z
        ECI z-location [meters]

    geodetic2eci() a.k.a lla2eci()
    """

    x, y, z = geodetic2ecef(lat, lon, alt, ell, deg)  # type: ignore[misc, arg-type]

    try:
        return ecef2eci(x, y, z, t)
    except NameError:
        raise ImportError("pip install numpy")


@overload
def enu2ecef(
    e1: float,
    n1: float,
    u1: float,
    lat0: float,
    lon0: float,
    h0: float,
    ell: Ellipsoid | None = None,
    deg: bool = True,
) -> tuple[float, float, float]:
    pass


@overload
def enu2ecef(
    e1: ArrayLike,
    n1: ArrayLike,
    u1: ArrayLike,
    lat0: ArrayLike,
    lon0: ArrayLike,
    h0: ArrayLike,
    ell: Ellipsoid | None = None,
    deg: bool = True,
) -> tuple[NDArray[Any], NDArray[Any], NDArray[Any]]:
    pass


def enu2ecef(
    e1: float | ArrayLike,
    n1: float | ArrayLike,
    u1: float | ArrayLike,
    lat0: float | ArrayLike,
    lon0: float | ArrayLike,
    h0: float | ArrayLike,
    ell: Ellipsoid | None = None,
    deg: bool = True,
) -> tuple[float, float, float] | tuple[NDArray[Any], NDArray[Any], NDArray[Any]]:
    """
    ENU to ECEF

    Parameters
    ----------

    e1
        target east ENU coordinate (meters)
    n1
        target north ENU coordinate (meters)
    u1
        target up ENU coordinate (meters)
    lat0
        Observer geodetic latitude
    lon0
        Observer geodetic longitude
    h0
         observer altitude above geodetic ellipsoid (meters)
    ell : Ellipsoid, optional
          reference ellipsoid
    deg : bool, optional
          degrees input/output  (False: radians in/out)


    Results
    -------
    x
        target x ECEF coordinate (meters)
    y
        target y ECEF coordinate (meters)
    z
        target z ECEF coordinate (meters)
    """
    x0, y0, z0 = geodetic2ecef(lat0, lon0, h0, ell, deg=deg)  # type: ignore[misc, arg-type]
    dx, dy, dz = enu2uvw(e1, n1, u1, lat0, lon0, deg=deg)  # type: ignore[misc, arg-type]

    return x0 + dx, y0 + dy, z0 + dz
