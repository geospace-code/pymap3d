"""Transforms involving ECEF: earth-centered, earth-fixed frame
as well as related ENU (east, north, up) and ECI (earth-centered inertial) frames.
"""

from __future__ import annotations

import warnings

try:
    from numpy import empty_like, finfo, where
except ImportError:
    pass

from datetime import datetime
from math import pi

from ._typing import FloatArray, FloatLike, NDArray
from .eci import ecef2eci, eci2ecef
from .ellipsoid import Ellipsoid, resolve_ellipsoid
from .mathfun import atan, atan2, cos, degrees, hypot, isclose, radians, sin, sqrt, tan

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


def geodetic2ecef(
    lat: FloatLike | FloatArray,
    lon: FloatLike | FloatArray,
    alt: FloatLike | FloatArray,
    ell: Ellipsoid | None = None,
    deg: bool = True,
) -> tuple:
    """
    point transformation from Geodetic of specified ellipsoid (default WGS-84) to ECEF

    Parameters
    ----------

    lat : array-like float
          point geodetic latitude
    lon : array-like float
          point geodetic longitude
    alt : array-like float
          point altitude above geodetic ellipsoid (meters)
    ell : Ellipsoid, optional
          reference ellipsoid
    deg : bool, optional
          degrees input/output  (False: radians in/out)


    Returns
    -------

    ECEF (Earth centered, Earth fixed)  x,y,z

    x : array-like float
        point x ECEF coordinate (meters)
    y : array-like float
        point y ECEF coordinate (meters)
    z : array-like float
        point z ECEF coordinate (meters)
    """

    if deg:
        lat = radians(lat)
        lon = radians(lon)

    ell = resolve_ellipsoid(ell)

    # radius of curvature of the prime vertical section
    N = ell.semimajor_axis**2 / hypot(ell.semimajor_axis * cos(lat), ell.semiminor_axis * sin(lat))
    # Compute cartesian (geocentric) coordinates given (curvilinear) geodetic coordinates.
    x = (N + alt) * cos(lat) * cos(lon)
    y = (N + alt) * cos(lat) * sin(lon)
    z = (N * (ell.semiminor_axis / ell.semimajor_axis) ** 2 + alt) * sin(lat)

    return x, y, z


def ecef2geodetic(
    x: FloatLike | FloatArray,
    y: FloatLike | FloatArray,
    z: FloatLike | FloatArray,
    ell: Ellipsoid | None = None,
    deg: bool = True,
) -> tuple:
    """
    convert ECEF (meters) to geodetic coordinates

    Parameters
    ----------
    x : array-like float
        point x ECEF coordinate (meters)
    y : array-like float
        point y ECEF coordinate (meters)
    z : array-like float
        point z ECEF coordinate (meters)
    ell : Ellipsoid, optional
          reference ellipsoid
    deg : bool, optional
          degrees input/output  (False: radians in/out)

    Returns
    -------
    lat : array-like float
          point geodetic latitude
    lon : array-like float
          point geodetic longitude
    alt : array-like float
          point altitude above geodetic ellipsoid (meters)

    based on:
    You, Rey-Jer. (2000). Transformation of Cartesian to Geodetic Coordinates without Iterations.
    Journal of Surveying Engineering. doi: 10.1061/(ASCE)0733-9453
    """

    ell = resolve_ellipsoid(ell)

    r = sqrt(x**2 + y**2 + z**2)

    E = sqrt(ell.semimajor_axis**2 - ell.semiminor_axis**2)

    # eqn. 4a
    u = sqrt(0.5 * (r**2 - E**2) + 0.5 * hypot(r**2 - E**2, 2 * E * z))

    hxy = hypot(x, y)

    huE = hypot(u, E)

    def _inside_numpy(x, y, z):
        # inside ellipsoid?
        return (
            x**2 / ell.semimajor_axis**2
            + y**2 / ell.semimajor_axis**2
            + z**2 / ell.semiminor_axis**2
            < 1
        )

    def _inside_scalar(x: float, y: float, z: float) -> bool:
        # inside ellipsoid?
        return (
            x**2 / ell.semimajor_axis**2
            + y**2 / ell.semimajor_axis**2
            + z**2 / ell.semiminor_axis**2
            < 1
        )

    def _lat_alt_numpy(x, y, z, Beta) -> tuple:
        """eqation 4c using Numpy"""
        lat = atan(ell.semimajor_axis / ell.semiminor_axis * tan(Beta))
        # patch latitude for float32 precision loss
        lim_pi2 = pi / 2 - finfo(Beta.dtype).eps
        lat = where(Beta >= lim_pi2, pi / 2, lat)
        lat = where(Beta <= -lim_pi2, -pi / 2, lat)

        # eqn. 7
        cosBeta = cos(Beta)
        # patch altitude for float32 precision loss
        cosBeta = where(Beta >= lim_pi2, 0, cosBeta)
        cosBeta = where(Beta <= -lim_pi2, 0, cosBeta)
        alt = hypot(z - ell.semiminor_axis * sin(Beta), hxy - ell.semimajor_axis * cosBeta)

        inside = _inside_numpy(x, y, z)
        try:
            alt[inside] = -alt[inside]
        except TypeError:
            if inside:
                alt = -alt

        return lat, alt

    def _beta_numpy(x, y, z) -> tuple:
        """using Numpy"""
        Beta = empty_like(r)
        ibad = isclose(u, 0) | isclose(hxy, 0)
        Beta[~ibad] = atan(huE[~ibad] / u[~ibad] * z[~ibad] / hxy[~ibad])
        # eqn. 13
        Beta[~ibad] += (
            (ell.semiminor_axis * u[~ibad] - ell.semimajor_axis * huE[~ibad] + E**2)
            * sin(Beta[~ibad])
        ) / (ell.semimajor_axis * huE[~ibad] * 1 / cos(Beta[~ibad]) - E**2 * cos(Beta[~ibad]))
        iz = ibad & isclose(z, 0)
        i1 = ibad & ~iz & (z > 0)
        i2 = ibad & ~iz & ~i1

        Beta[iz] = 0
        Beta[i1] = pi / 2
        Beta[i2] = -pi / 2

        lat, alt = _lat_alt_numpy(x, y, z, Beta)

        lat = lat.squeeze()[()]

        return lat, alt

    def _beta_scalar(x: float, y: float, z: float) -> tuple:
        """using scalar math"""
        try:
            with warnings.catch_warnings(record=True):
                warnings.simplefilter("error")
                Beta = atan(huE / u * z / hxy)
                # eqn. 13
                Beta += (
                    (ell.semiminor_axis * u - ell.semimajor_axis * huE + E**2) * sin(Beta)
                ) / (ell.semimajor_axis * huE * 1 / cos(Beta) - E**2 * cos(Beta))
        except (ArithmeticError, RuntimeWarning):
            if isclose(z, 0):
                Beta = 0
            elif z > 0:
                Beta = pi / 2
            else:
                Beta = -pi / 2

        # eqn. 4c
        lat = atan(ell.semimajor_axis / ell.semiminor_axis * tan(Beta))
        # eqn. 7
        cosBeta = cos(Beta)

        alt = hypot(z - ell.semiminor_axis * sin(Beta), hxy - ell.semimajor_axis * cosBeta)

        if _inside_scalar(x, y, z):
            alt = -alt

        return lat, alt

    # eqn. 4b
    try:
        lat, alt = _beta_numpy(x, y, z)
    except (NameError, TypeError):
        lat, alt = _beta_scalar(float(x), float(y), float(z))

    lon = atan2(y, x)

    if deg:
        lat = degrees(lat)
        lon = degrees(lon)

    return lat, lon, alt


def ecef2enuv(
    u: FloatLike, v: FloatLike, w: FloatLike, lat0: FloatLike, lon0: FloatLike, deg: bool = True
) -> tuple:
    """
    VECTOR from observer to target  ECEF => ENU

    Parameters
    ----------
    u : array-like float
        target x ECEF coordinate (meters)
    v : array-like float
        target y ECEF coordinate (meters)
    w : array-like float
        target z ECEF coordinate (meters)
    lat0 : array-like float
        Observer geodetic latitude
    lon0 : array-like float
        Observer geodetic longitude
    deg : bool, optional
        degrees input/output  (False: radians in/out)

    Returns
    -------
    uEast : array-like float
        target east ENU coordinate (meters)
    vNorth : array-like float
        target north ENU coordinate (meters)
    wUp : array-like float
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


def ecef2enu(
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
    from observer to target, ECEF => ENU

    Parameters
    ----------
    x : array-like float
        target x ECEF coordinate (meters)
    y : array-like float
        target y ECEF coordinate (meters)
    z : array-like float
        target z ECEF coordinate (meters)
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

    Returns
    -------
    East : array-like float
        target east ENU coordinate (meters)
    North : array-like float
        target north ENU coordinate (meters)
    Up: array-like float
        target up ENU coordinate (meters)
    """
    x0, y0, z0 = geodetic2ecef(lat0, lon0, h0, ell, deg=deg)

    return uvw2enu(x - x0, y - y0, z - z0, lat0, lon0, deg=deg)


def enu2uvw(
    east: FloatLike,
    north: FloatLike,
    up: FloatLike,
    lat0: FloatLike,
    lon0: FloatLike,
    deg: bool = True,
) -> tuple:
    """
    enu2uvw gives the relative ECEF displacement vector based on ENU position for rotation only.
    Common in radar, interferometry, baseline calculations, etc.

    See enu2ecef() for observer to target, which includes rotation and translation.

    Parameters
    ----------

    east : array-like float
        target east ENU coordinate (meters)
    north : array-like float
        target north ENU coordinate (meters)
    up : array-like float
        target up ENU coordinate (meters)
    lat0 : array-like float
        Observer geodetic latitude
    lon0 : array-like float
        Observer geodetic longitude
    deg : bool, optional
        degrees input/output  (False: radians in/out)

    Results
    -------

    u : array-like float
        x ECEF displacement (meters)
    v : array-like float
        y ECEF displacement (meters)
    w : array-like float
        z ECEF displacement (meters)
    """

    if deg:
        lat0 = radians(lat0)
        lon0 = radians(lon0)

    t = cos(lat0) * up - sin(lat0) * north
    w = sin(lat0) * up + cos(lat0) * north

    u = cos(lon0) * t - sin(lon0) * east
    v = sin(lon0) * t + cos(lon0) * east

    return u, v, w


def uvw2enu(
    u: FloatLike, v: FloatLike, w: FloatLike, lat0: FloatLike, lon0: FloatLike, deg: bool = True
) -> tuple:
    """
    uvw2enu gives the relative ENU position based on ECEF displacement vector for rotation only.
    Common in radar, interferometry, baseline calculations, etc.

    See ecef2enu() for observer to target, which includes rotation and translation.

    Parameters
    ----------

    u : array-like float
        target x ECEF coordinate (meters)
    v : array-like float
        target y ECEF coordinate (meters)
    w : array-like float
        target z ECEF coordinate (meters)
    lat0 : array-like float
        Observer geodetic latitude
    lon0 : array-like float
        Observer geodetic longitude
    deg : bool, optional
        degrees input/output  (False: radians in/out)

    Results
    -------

    East : array-like float
        East ENU coordinate (meters)
    North : array-like float
        North ENU coordinate (meters)
    Up : array-like float
        Up ENU coordinate (meters)
    """
    if deg:
        lat0 = radians(lat0)
        lon0 = radians(lon0)

    t = cos(lon0) * u + sin(lon0) * v
    East = -sin(lon0) * u + cos(lon0) * v
    Up = cos(lat0) * t + sin(lat0) * w
    North = -sin(lat0) * t + cos(lat0) * w

    return East, North, Up


def eci2geodetic(
    x: FloatLike,
    y: FloatLike,
    z: FloatLike,
    t: datetime,
    ell: Ellipsoid | None = None,
    *,
    deg: bool = True,
) -> tuple:
    """
    convert Earth Centered Internal ECI to geodetic coordinates

    J2000 time

    Parameters
    ----------
    x : array-like float
        ECI x-location [meters]
    y : array-like float
        ECI y-location [meters]
    z : array-like float
        ECI z-location [meters]
    t : datetime.datetime, float
        UTC time
    ell : Ellipsoid, optional
        planet ellipsoid model
    deg : bool, optional
        if True, degrees. if False, radians

    Results
    -------
    lat : array-like float
        geodetic latitude
    lon : array-like float
        geodetic longitude
    alt : array-like float
        altitude above ellipsoid  (meters)

    eci2geodetic() is also known as eci2lla() by other packages e.g. MATLAB Aerospace Toolbox.
    """

    xecef, yecef, zecef = eci2ecef(x, y, z, t)

    return ecef2geodetic(xecef, yecef, zecef, ell, deg)


def geodetic2eci(
    lat: FloatLike,
    lon: FloatLike,
    alt: FloatLike,
    t: datetime,
    ell: Ellipsoid | None = None,
    *,
    deg: bool = True,
) -> tuple:
    """
    convert geodetic coordinates to Earth Centered Internal ECI

    J2000 frame

    Parameters
    ----------
    lat : array-like float
        geodetic latitude
    lon : array-like float
        geodetic longitude
    alt : array-like float
        altitude above ellipsoid  (meters)
    t : datetime.datetime, float
        UTC time
    ell : Ellipsoid, optional
        planet ellipsoid model
    deg : bool, optional
        if True, degrees. if False, radians

    Results
    -------
    x : array-like float
        ECI x-location [meters]
    y : array-like float
        ECI y-location [meters]
    z : array-like float
        ECI z-location [meters]

    geodetic2eci() a.k.a lla2eci()
    """

    x, y, z = geodetic2ecef(lat, lon, alt, ell, deg)

    return ecef2eci(x, y, z, t)


def enu2ecef(
    e1: FloatLike,
    n1: FloatLike,
    u1: FloatLike,
    lat0: FloatLike,
    lon0: FloatLike,
    h0: FloatLike,
    ell: Ellipsoid | None = None,
    deg: bool = True,
) -> tuple:
    """
    ENU to ECEF

    Parameters
    ----------

    e1 : array-like float
        target east ENU coordinate (meters)
    n1 : array-like float
        target north ENU coordinate (meters)
    u1 : array-like float
        target up ENU coordinate (meters)
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
    x : array-like float
        target x ECEF coordinate (meters)
    y : array-like float
        target y ECEF coordinate (meters)
    z : array-like float
        target z ECEF coordinate (meters)
    """

    x0, y0, z0 = geodetic2ecef(lat0, lon0, h0, ell, deg=deg)
    u, v, w = enu2uvw(e1, n1, u1, lat0, lon0, deg=deg)

    return x0 + u, y0 + v, z0 + w
