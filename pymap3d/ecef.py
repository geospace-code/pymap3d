""" Transforms involving ECEF: earth-centered, earth-fixed frame """
try:
    from numpy import radians, sin, cos, tan, arctan as atan, hypot, degrees, arctan2 as atan2, sqrt, pi, vectorize
    from .eci import eci2ecef
except ImportError:
    from math import radians, sin, cos, tan, atan, hypot, degrees, atan2, sqrt, pi

    vectorize = eci2ecef = None
import typing
from datetime import datetime

from .ellipsoid import Ellipsoid
from .utils import sanitize

# py < 3.6 compatible
tau = 2 * pi

__all__ = ["geodetic2ecef", "ecef2geodetic", "ecef2enuv", "ecef2enu", "enu2uvw", "uvw2enu", "eci2geodetic", "enu2ecef"]


def geodetic2ecef(lat: float, lon: float, alt: float, ell: Ellipsoid = None, deg: bool = True) -> typing.Tuple[float, float, float]:
    """
    point transformation from Geodetic of specified ellipsoid (default WGS-84) to ECEF

    Parameters
    ----------

    lat : float
           target geodetic latitude
    lon : float
           target geodetic longitude
    h : float
         target altitude above geodetic ellipsoid (meters)
    ell : Ellipsoid, optional
          reference ellipsoid
    deg : bool, optional
          degrees input/output  (False: radians in/out)


    Returns
    -------

    ECEF (Earth centered, Earth fixed)  x,y,z

    x : float
        target x ECEF coordinate (meters)
    y : float
        target y ECEF coordinate (meters)
    z : float
        target z ECEF coordinate (meters)
    """
    lat, ell = sanitize(lat, ell, deg)
    if deg:
        lon = radians(lon)

    # radius of curvature of the prime vertical section
    N = ell.semimajor_axis ** 2 / sqrt(ell.semimajor_axis ** 2 * cos(lat) ** 2 + ell.semiminor_axis ** 2 * sin(lat) ** 2)
    # Compute cartesian (geocentric) coordinates given  (curvilinear) geodetic
    # coordinates.
    x = (N + alt) * cos(lat) * cos(lon)
    y = (N + alt) * cos(lat) * sin(lon)
    z = (N * (ell.semiminor_axis / ell.semimajor_axis) ** 2 + alt) * sin(lat)

    return x, y, z


def ecef2geodetic(x: float, y: float, z: float, ell: Ellipsoid = None, deg: bool = True) -> typing.Tuple[float, float, float]:
    if vectorize is not None:
        fun = vectorize(ecef2geodetic_point)
        return fun(x, y, z, ell, deg)
    else:
        return ecef2geodetic_point(x, y, z, ell, deg)


def ecef2geodetic_point(x: float, y: float, z: float, ell: Ellipsoid = None, deg: bool = True) -> typing.Tuple[float, float, float]:
    """
    convert ECEF (meters) to geodetic coordinates

    Parameters
    ----------
    x : float
        target x ECEF coordinate (meters)
    y : float
        target y ECEF coordinate (meters)
    z : float
        target z ECEF coordinate (meters)
    ell : Ellipsoid, optional
          reference ellipsoid
    deg : bool, optional
          degrees input/output  (False: radians in/out)

    Returns
    -------
    lat : float
           target geodetic latitude
    lon : float
           target geodetic longitude
    h : float
         target altitude above geodetic ellipsoid (meters)

    based on:
    You, Rey-Jer. (2000). Transformation of Cartesian to Geodetic Coordinates without Iterations.
    Journal of Surveying Engineering. doi: 10.1061/(ASCE)0733-9453
    """
    if ell is None:
        ell = Ellipsoid()

    r = sqrt(x ** 2 + y ** 2 + z ** 2)

    E = sqrt(ell.semimajor_axis ** 2 - ell.semiminor_axis ** 2)

    # eqn. 4a
    u = sqrt(0.5 * (r ** 2 - E ** 2) + 0.5 * sqrt((r ** 2 - E ** 2) ** 2 + 4 * E ** 2 * z ** 2))

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
    eps = ((ell.semiminor_axis * u - ell.semimajor_axis * huE + E ** 2) * sin(Beta)) / (
        ell.semimajor_axis * huE * 1 / cos(Beta) - E ** 2 * cos(Beta)
    )

    Beta += eps
    # %% final output
    lat = atan(ell.semimajor_axis / ell.semiminor_axis * tan(Beta))

    lon = atan2(y, x)

    # eqn. 7
    alt = hypot(z - ell.semiminor_axis * sin(Beta), Q - ell.semimajor_axis * cos(Beta))

    # inside ellipsoid?
    inside = x ** 2 / ell.semimajor_axis ** 2 + y ** 2 / ell.semimajor_axis ** 2 + z ** 2 / ell.semiminor_axis ** 2 < 1
    if inside:
        alt = -alt

    if deg:
        lat = degrees(lat)
        lon = degrees(lon)

    return lat, lon, alt


def ecef2enuv(u: float, v: float, w: float, lat0: float, lon0: float, deg: bool = True) -> typing.Tuple[float, float, float]:
    """
    VECTOR from observer to target  ECEF => ENU

    Parameters
    ----------
    u : float
        target x ECEF coordinate (meters)
    v : float
        target y ECEF coordinate (meters)
    w : float
        target z ECEF coordinate (meters)
    lat0 : float
           Observer geodetic latitude
    lon0 : float
           Observer geodetic longitude
    h0 : float
         observer altitude above geodetic ellipsoid (meters)
    deg : bool, optional
          degrees input/output  (False: radians in/out)

    Returns
    -------
    uEast : float
        target east ENU coordinate (meters)
    vNorth : float
        target north ENU coordinate (meters)
    wUp : float
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
    x: float, y: float, z: float, lat0: float, lon0: float, h0: float, ell: Ellipsoid = None, deg: bool = True
) -> typing.Tuple[float, float, float]:
    """
    from observer to target, ECEF => ENU

    Parameters
    ----------
    x : float
        target x ECEF coordinate (meters)
    y : float
        target y ECEF coordinate (meters)
    z : float
        target z ECEF coordinate (meters)
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
    East : float
        target east ENU coordinate (meters)
    North : float
        target north ENU coordinate (meters)
    Up : float
        target up ENU coordinate (meters)

    """
    x0, y0, z0 = geodetic2ecef(lat0, lon0, h0, ell, deg=deg)

    return uvw2enu(x - x0, y - y0, z - z0, lat0, lon0, deg=deg)


def enu2uvw(east: float, north: float, up: float, lat0: float, lon0: float, deg: bool = True) -> typing.Tuple[float, float, float]:
    """
    Parameters
    ----------

    e1 : float
        target east ENU coordinate (meters)
    n1 : float
        target north ENU coordinate (meters)
    u1 : float
        target up ENU coordinate (meters)

    Results
    -------

    u : float
    v : float
    w : float
    """

    if deg:
        lat0 = radians(lat0)
        lon0 = radians(lon0)

    t = cos(lat0) * up - sin(lat0) * north
    w = sin(lat0) * up + cos(lat0) * north

    u = cos(lon0) * t - sin(lon0) * east
    v = sin(lon0) * t + cos(lon0) * east

    return u, v, w


def uvw2enu(u: float, v: float, w: float, lat0: float, lon0: float, deg: bool = True) -> typing.Tuple[float, float, float]:
    """
    Parameters
    ----------

    u : float
    v : float
    w : float


    Results
    -------

    East : float
        target east ENU coordinate (meters)
    North : float
        target north ENU coordinate (meters)
    Up : float
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


def eci2geodetic(x: float, y: float, z: float, t: datetime, useastropy: bool = True) -> typing.Tuple[float, float, float]:
    """
    convert ECI to geodetic coordinates

    Parameters
    ----------
    x : float
        ECI x-location [meters]
    y : float
        ECI y-location [meters]
    z : float
        ECI z-location [meters]
    t : datetime.datetime, float
          length N vector of datetime OR greenwich sidereal time angle [radians].

    Results
    -------
    lat : float
          geodetic latitude
    lon : float
          geodetic longitude
    alt : float
          altitude above ellipsoid  (meters)

    Notes
    -----

    Conversion is idealized: doesn't consider nutations, perterbations,
    etc. like the IAU-76/FK5 or IAU-2000/2006 model-based conversions
    from ECI to ECEF

    eci2geodetic() a.k.a. eci2lla()
    """
    if eci2ecef is None:
        raise ImportError("pip install numpy")

    xecef, yecef, zecef = eci2ecef(x, y, z, t, useastropy=useastropy)

    return ecef2geodetic(xecef, yecef, zecef)


def enu2ecef(
    e1: float, n1: float, u1: float, lat0: float, lon0: float, h0: float, ell: Ellipsoid = None, deg: bool = True
) -> typing.Tuple[float, float, float]:
    """
    ENU to ECEF

    Parameters
    ----------

    e1 : float
        target east ENU coordinate (meters)
    n1 : float
        target north ENU coordinate (meters)
    u1 : float
        target up ENU coordinate (meters)
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
        target x ECEF coordinate (meters)
    y : float
        target y ECEF coordinate (meters)
    z : float
        target z ECEF coordinate (meters)
    """
    x0, y0, z0 = geodetic2ecef(lat0, lon0, h0, ell, deg=deg)
    dx, dy, dz = enu2uvw(e1, n1, u1, lat0, lon0, deg=deg)

    return x0 + dx, y0 + dy, z0 + dz
