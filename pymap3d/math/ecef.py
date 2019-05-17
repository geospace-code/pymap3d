""" Transforms involving ECEF: earth-centered, earth-fixed frame """
from math import radians, sin, cos, tan, hypot, degrees, sqrt, pi
from math import atan as arctan, atan2 as arctan2
from ..ellipsoid import Ellipsoid

tau = 2 * pi

__all__ = ['geodetic2ecef', 'ecef2geodetic', 'ecef2enuv', 'ecef2enu', 'enu2uvw', 'uvw2enu', 'enu2ecef']


# def geodetic2ecef(lat: float, lon: float, alt: float, ell: Ellipsoid = None, deg: bool = True):
def geodetic2ecef(lat, lon, alt, ell=None, deg=True):
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
    if ell is None:
        ell = Ellipsoid()

    if deg:
        lat = radians(lat)
        lon = radians(lon)

    if (lat < -pi / 2) | (lat > pi / 2):
        raise ValueError('-90 <= lat <= 90')

    # radius of curvature of the prime vertical section
    N = ell.semimajor_axis**2 / sqrt(ell.semimajor_axis**2 * cos(lat)**2 + ell.semiminor_axis**2 * sin(lat)**2)
    # Compute cartesian (geocentric) coordinates given  (curvilinear) geodetic
    # coordinates.
    x = (N + alt) * cos(lat) * cos(lon)
    y = (N + alt) * cos(lat) * sin(lon)
    z = (N * (ell.semiminor_axis / ell.semimajor_axis)**2 + alt) * sin(lat)

    return x, y, z


# def ecef2geodetic(x: float, y: float, z: float, ell: Ellipsoid = None, deg: bool = True):
def ecef2geodetic(x, y, z, ell=None, deg=True):
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

    r = sqrt(x**2 + y**2 + z**2)

    E = sqrt(ell.semimajor_axis**2 - ell.semiminor_axis**2)

    # eqn. 4a
    u = sqrt(0.5 * (r**2 - E**2) + 0.5 * sqrt((r**2 - E**2)**2 + 4 * E**2 * z**2))

    Q = hypot(x, y)

    huE = hypot(u, E)

    # eqn. 4b
    Beta = arctan(huE / u * z / hypot(x, y))

    # eqn. 13
    eps = ((ell.semiminor_axis * u - ell.semimajor_axis * huE + E**2) * sin(Beta)) / (ell.semimajor_axis * huE * 1 / cos(Beta) - E**2 * cos(Beta))

    Beta += eps
# %% final output
    lat = arctan(ell.semimajor_axis / ell.semiminor_axis * tan(Beta))

    lon = arctan2(y, x)

    # eqn. 7
    alt = hypot(z - ell.semiminor_axis * sin(Beta),
                Q - ell.semimajor_axis * cos(Beta))

    # inside ellipsoid?
    inside = x**2 / ell.semimajor_axis**2 + y**2 / ell.semimajor_axis**2 + z**2 / ell.semiminor_axis**2 < 1
    if inside:
        alt = -alt

    if deg:
        lat = degrees(lat)
        lon = degrees(lon)

    return lat, lon, alt


# def ecef2enuv(u: float, v: float, w: float, lat0: float, lon0: float, deg: bool = True):
def ecef2enuv(u, v, w, lat0, lon0, deg=True):
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


# def ecef2enu(x: float, y: float, z: float, lat0: float, lon0: float, h0: float, ell: Ellipsoid = None, deg: bool = True):
def ecef2enu(x, y, z, lat0, lon0, h0, ell=None, deg=True):
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


# def enu2uvw(east: float, north: float, up: float, lat0: float, lon0: float, deg: bool = True):
def enu2uvw(east, north, up, lat0, lon0, deg=True):
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


# def uvw2enu(u: float, v: float, w: float, lat0: float, lon0: float, deg: bool = True):
def uvw2enu(u, v, w, lat0, lon0, deg=True):
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


# def enu2ecef(e1: float, n1: float, u1: float, lat0: float, lon0: float, h0: float, ell: Ellipsoid = None, deg: bool = True):
def enu2ecef(e1, n1, u1, lat0, lon0, h0, ell=None, deg=True):
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
