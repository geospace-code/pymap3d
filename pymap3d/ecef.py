""" Transforms involving ECEF: earth-centered, earth-fixed frame """
from numpy import radians, sin, cos, tan, arctan, hypot, degrees, arctan2, sqrt, pi
import numpy as np
from typing import Tuple
from datetime import datetime
try:
    from math import tau  # py >= 3.6
except ImportError:
    tau = 2 * pi

from .eci import eci2ecef

__all__ = ['Ellipsoid', 'geodetic2ecef', 'ecef2geodetic', 'ecef2enuv', 'ecef2enu', 'enu2uvw', 'uvw2enu', 'eci2geodetic', 'enu2ecef']


class Ellipsoid:
    """
    generate reference ellipsoid parameters

    https://en.wikibooks.org/wiki/PROJ.4#Spheroid

    https://tharsis.gsfc.nasa.gov/geodesy.html

    https://nssdc.gsfc.nasa.gov/planetary/factsheet/index.html
    """

    def __init__(self, model: str = 'wgs84'):
        """
        feel free to suggest additional ellipsoids

        Parameters
        ----------
        model : str
                name of ellipsoid
        """
        if model == 'wgs84':
            """https://en.wikipedia.org/wiki/World_Geodetic_System#WGS84"""
            self.a = 6378137.  # semi-major axis [m]
            self.f = 1 / 298.2572235630  # flattening
            self.b = self.a * (1 - self.f)  # semi-minor axis
        elif model == 'grs80':
            """https://en.wikipedia.org/wiki/GRS_80"""
            self.a = 6378137.  # semi-major axis [m]
            self.f = 1 / 298.257222100882711243  # flattening
            self.b = self.a * (1 - self.f)  # semi-minor axis
        elif model == 'clrk66':  # Clarke 1866
            self.a = 6378206.4  # semi-major axis [m]
            self.b = 6356583.8  # semi-minor axis
            self.f = -(self.b / self.a - 1)
        elif model == 'mars':  # https://tharsis.gsfc.nasa.gov/geodesy.html
            self.a = 3396900
            self.b = 3376097.80585952
            self.f = 1 / 163.295274386012
        elif model == 'moon':
            self.a = 1738000.
            self.b = self.a
            self.f = 0.
        elif model == 'venus':
            self.a = 6051000.
            self.b = self.a
            self.f = 0.
        elif model == 'pluto':
            self.a = 1187000.
            self.b = self.a
            self.f = 0.
        else:
            raise NotImplementedError('{} model not implemented, let us know and we will add it (or make a pull request)'.format(model))


def get_radius_normal(lat_radians: float, ell: Ellipsoid = None) -> float:
    """
    Compute normal radius of planetary body

    Parameters
    ----------

    lat_radians : float
        latitude in radians
    ell : Ellipsoid, optional
        reference ellipsoid

    Returns
    -------

    radius : float
        normal radius (meters)
    """
    if ell is None:
        ell = Ellipsoid()

    a = ell.a
    b = ell.b

    return a**2 / sqrt(a**2 * cos(lat_radians)**2 + b**2 * sin(lat_radians)**2)


def geodetic2ecef(lat: float, lon: float, alt: float,
                  ell: Ellipsoid = None, deg: bool = True) -> Tuple[float, float, float]:
    """
    point transformation from Geodetic of specified ellipsoid (default WGS-84) to ECEF

    Parameters
    ----------

    lat : float or numpy.ndarray of float
           target geodetic latitude
    lon : float or numpy.ndarray of float
           target geodetic longitude
    h : float or numpy.ndarray of float
         target altitude above geodetic ellipsoid (meters)
    ell : Ellipsoid, optional
          reference ellipsoid
    deg : bool, optional
          degrees input/output  (False: radians in/out)


    Returns
    -------

    ECEF (Earth centered, Earth fixed)  x,y,z

    x : float or numpy.ndarray of float
        target x ECEF coordinate (meters)
    y : float or numpy.ndarray of float
        target y ECEF coordinate (meters)
    z : float or numpy.ndarray of float
        target z ECEF coordinate (meters)
    """
    if ell is None:
        ell = Ellipsoid()

    if deg:
        lat = radians(lat)
        lon = radians(lon)

    with np.errstate(invalid='ignore'):
        # need np.any() to handle scalar and array cases
        if np.any((lat < -pi / 2) | (lat > pi / 2)):
            raise ValueError('-90 <= lat <= 90')

    # radius of curvature of the prime vertical section
    N = get_radius_normal(lat, ell)
    # Compute cartesian (geocentric) coordinates given  (curvilinear) geodetic
    # coordinates.
    x = (N + alt) * cos(lat) * cos(lon)
    y = (N + alt) * cos(lat) * sin(lon)
    z = (N * (ell.b / ell.a)**2 + alt) * sin(lat)

    return x, y, z


def ecef2geodetic(x: float, y: float, z: float,
                  ell: Ellipsoid = None, deg: bool = True) -> Tuple[float, float, float]:
    """
    convert ECEF (meters) to geodetic coordinates

    Parameters
    ----------
    x : float or numpy.ndarray of float
        target x ECEF coordinate (meters)
    y : float or numpy.ndarray of float
        target y ECEF coordinate (meters)
    z : float or numpy.ndarray of float
        target z ECEF coordinate (meters)
    ell : Ellipsoid, optional
          reference ellipsoid
    deg : bool, optional
          degrees input/output  (False: radians in/out)

    Returns
    -------
    lat : float or numpy.ndarray of float
           target geodetic latitude
    lon : float or numpy.ndarray of float
           target geodetic longitude
    h : float or numpy.ndarray of float
         target altitude above geodetic ellipsoid (meters)

    based on:
    You, Rey-Jer. (2000). Transformation of Cartesian to Geodetic Coordinates without Iterations.
    Journal of Surveying Engineering. doi: 10.1061/(ASCE)0733-9453
    """
    if ell is None:
        ell = Ellipsoid()

    x = np.asarray(x)
    y = np.asarray(y)
    z = np.asarray(z)

    r = sqrt(x**2 + y**2 + z**2)

    E = sqrt(ell.a**2 - ell.b**2)

    # eqn. 4a
    u = sqrt(0.5 * (r**2 - E**2) + 0.5 * sqrt((r**2 - E**2)**2 + 4 * E**2 * z**2))

    Q = hypot(x, y)

    huE = hypot(u, E)

    # eqn. 4b
    with np.errstate(divide='ignore'):
        Beta = arctan(huE / u * z / hypot(x, y))

    # eqn. 13
    eps = ((ell.b * u - ell.a * huE + E**2) * sin(Beta)) / (ell.a * huE * 1 / cos(Beta) - E**2 * cos(Beta))

    Beta += eps
# %% final output
    lat = arctan(ell.a / ell.b * tan(Beta))

    lon = arctan2(y, x)

    # eqn. 7
    alt = hypot(z - ell.b * sin(Beta),
                Q - ell.a * cos(Beta))

    # inside ellipsoid?
    with np.errstate(invalid='ignore'):
        inside = x**2 / ell.a**2 + y**2 / ell.a**2 + z**2 / ell.b**2 < 1
    if isinstance(inside, np.ndarray):
        alt[inside] = -alt[inside]
    elif inside:
        alt = -alt

    if deg:
        lat = degrees(lat)
        lon = degrees(lon)

    return lat, lon, alt


def ecef2enuv(u: float, v: float, w: float,
              lat0: float, lon0: float, deg: bool = True) -> Tuple[float, float, float]:
    """
    VECTOR from observer to target  ECEF => ENU

    Parameters
    ----------
    u : float or numpy.ndarray of float
        target x ECEF coordinate (meters)
    v : float or numpy.ndarray of float
        target y ECEF coordinate (meters)
    w : float or numpy.ndarray of float
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
    uEast : float or numpy.ndarray of float
        target east ENU coordinate (meters)
    vNorth : float or numpy.ndarray of float
        target north ENU coordinate (meters)
    wUp : float or numpy.ndarray of float
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


def ecef2enu(x: float, y: float, z: float,
             lat0: float, lon0: float, h0: float,
             ell: Ellipsoid = None, deg: bool = True) -> Tuple[float, float, float]:
    """
    from observer to target, ECEF => ENU

    Parameters
    ----------
    x : float or numpy.ndarray of float
        target x ECEF coordinate (meters)
    y : float or numpy.ndarray of float
        target y ECEF coordinate (meters)
    z : float or numpy.ndarray of float
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
    East : float or numpy.ndarray of float
        target east ENU coordinate (meters)
    North : float or numpy.ndarray of float
        target north ENU coordinate (meters)
    Up : float or numpy.ndarray of float
        target up ENU coordinate (meters)

    """
    x0, y0, z0 = geodetic2ecef(lat0, lon0, h0, ell, deg=deg)

    return uvw2enu(x - x0, y - y0, z - z0, lat0, lon0, deg=deg)


def enu2uvw(east: float, north: float, up: float,
            lat0: float, lon0: float, deg: bool = True) -> Tuple[float, float, float]:
    """
    Parameters
    ----------

    e1 : float or numpy.ndarray of float
        target east ENU coordinate (meters)
    n1 : float or numpy.ndarray of float
        target north ENU coordinate (meters)
    u1 : float or numpy.ndarray of float
        target up ENU coordinate (meters)

    Results
    -------

    u : float or numpy.ndarray of float
    v : float or numpy.ndarray of float
    w : float or numpy.ndarray of float
    """

    if deg:
        lat0 = radians(lat0)
        lon0 = radians(lon0)

    t = cos(lat0) * up - sin(lat0) * north
    w = sin(lat0) * up + cos(lat0) * north

    u = cos(lon0) * t - sin(lon0) * east
    v = sin(lon0) * t + cos(lon0) * east

    return u, v, w


def uvw2enu(u: float, v: float, w: float,
            lat0: float, lon0: float, deg: bool = True) -> Tuple[float, float, float]:
    """
    Parameters
    ----------

    u : float or numpy.ndarray of float
    v : float or numpy.ndarray of float
    w : float or numpy.ndarray of float


    Results
    -------

    East : float or numpy.ndarray of float
        target east ENU coordinate (meters)
    North : float or numpy.ndarray of float
        target north ENU coordinate (meters)
    Up : float or numpy.ndarray of float
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


def eci2geodetic(eci: np.ndarray, t: datetime,
                 useastropy: bool = True) -> Tuple[float, float, float]:
    """
    convert ECI to geodetic coordinates

    Parameters
    ----------
    eci : tuple of float
          [meters] Nx3 target ECI location (x,y,z)
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
    ecef = np.atleast_2d(eci2ecef(eci, t, useastropy=useastropy))

    return np.asarray(ecef2geodetic(ecef[:, 0], ecef[:, 1], ecef[:, 2])).squeeze()


def enu2ecef(e1: float, n1: float, u1: float,
             lat0: float, lon0: float, h0: float,
             ell: Ellipsoid = None, deg: bool = True) -> Tuple[float, float, float]:
    """
    ENU to ECEF

    Parameters
    ----------

    e1 : float or numpy.ndarray of float
        target east ENU coordinate (meters)
    n1 : float or numpy.ndarray of float
        target north ENU coordinate (meters)
    u1 : float or numpy.ndarray of float
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
    x : float or numpy.ndarray of float
        target x ECEF coordinate (meters)
    y : float or numpy.ndarray of float
        target y ECEF coordinate (meters)
    z : float or numpy.ndarray of float
        target z ECEF coordinate (meters)
    """
    x0, y0, z0 = geodetic2ecef(lat0, lon0, h0, ell, deg=deg)
    dx, dy, dz = enu2uvw(e1, n1, u1, lat0, lon0, deg=deg)

    return x0 + dx, y0 + dy, z0 + dz
