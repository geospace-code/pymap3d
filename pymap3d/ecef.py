from numpy import radians, sin, cos, tan, allclose, hypot, degrees, arctan2, sqrt, pi
import numpy as np
from copy import deepcopy
from typing import Tuple
from datetime import datetime

from .eci import eci2ecef


class Ellipsoid:
    """
    generate reference ellipsoid parameters

    https://en.wikibooks.org/wiki/PROJ.4#Spheroid
    """

    def __init__(self, model: str='wgs84') -> None:
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
            self.b = 1738000.
            self.f = 0.
        elif model == 'venus':
            self.a = 6051000.
            self.b = 6051000.
            self.f = 0.
        else:
            raise NotImplementedError('{} model not implemented, let us know and we will add it (or make a pull request)'.format(model))


def get_radius_normal(lat_radians: float, ell: Ellipsoid=None) -> float:
    """ Compute normal radius of planetary body"""
    if ell is None:
        ell = Ellipsoid()

    a = ell.a
    b = ell.b

    return a**2 / sqrt(a**2 * cos(lat_radians)**2 + b**2 * sin(lat_radians)**2)


def geodetic2ecef(lat: float, lon: float, alt: float,
                  ell: Ellipsoid=None, deg: bool=True) -> Tuple[float, float, float]:
    """
    Point

    input:
    -----
    lat, lon (degrees)
    alt (altitude, meters)    [0, Infinity)
    ell    reference ellipsoid
    deg    degrees input/output  (False: radians in/out)

    output: ECEF x,y,z (meters)
    """
    if ell is None:
        ell = Ellipsoid()

    if deg:
        lat = radians(lat)
        lon = radians(lon)

    with np.errstate(invalid='ignore'):
        if ((lat < -pi / 2) | (lat > pi / 2)).any():
            raise ValueError('-90 <= lat <= 90')

        if ((lon < -pi) | (lon > 2 * pi)).any():
            raise ValueError('-180 <= lat <= 360')

        if (np.asarray(alt) < 0).any():
            raise ValueError('altitude \in  [0, Infinity)')
    # radius of curvature of the prime vertical section
    N = get_radius_normal(lat, ell)
    # Compute cartesian (geocentric) coordinates given  (curvilinear) geodetic
    # coordinates.
    x = (N + alt) * cos(lat) * cos(lon)
    y = (N + alt) * cos(lat) * sin(lon)
    z = (N * (ell.b / ell.a)**2 + alt) * sin(lat)

    return x, y, z


def ecef2geodetic(x: float, y: float, z: float,
                  ell: Ellipsoid=None, deg: bool=True) -> Tuple[float, float, float]:
    """
    convert ECEF (meters) to geodetic coordinates

    input
    -----
    x,y,z  [meters] target ECEF location                             [0,Infinity)
    ell    reference ellipsoid
    deg    degrees input/output  (False: radians in/out)

    output
    ------
    lat,lon   (degrees/radians)
    alt  (meters)

    Algorithm is based on
    http://www.astro.uni.torun.pl/~kb/Papers/geod/Geod-BG.htm
    This algorithm provides a converging solution to the latitude equation
    in terms of the parametric or reduced latitude form (v)
    This algorithm provides a uniform solution over all latitudes as it does
    not involve division by cos(phi) or sin(phi)
    """
    if ell is None:
        ell = Ellipsoid()

    ea = ell.a
    eb = ell.b
    rad = hypot(x, y)
# Constant required for Latitude equation
    rho = arctan2(eb * z, ea * rad)
# Constant required for latitude equation
    c = (ea**2 - eb**2) / hypot(ea * rad, eb * z)
# Starter for the Newtons Iteration Method
    vnew = arctan2(ea * z, eb * rad)
# Initializing the parametric latitude
    v = 0
    for _ in range(5):
        v = deepcopy(vnew)
# %% Newtons Method for computing iterations
        vnew = v - ((2 * sin(v - rho) - c * sin(2 * v)) /
                    (2 * (cos(v - rho) - c * cos(2 * v))))

        if allclose(v, vnew):
            break
# %% Computing latitude from the root of the latitude equation
    lat = arctan2(ea * tan(vnew), eb)
    # by inspection
    lon = arctan2(y, x)

    alt = (((rad - ea * cos(vnew)) * cos(lat)) +
           ((z - eb * sin(vnew)) * sin(lat)))

    with np.errstate(invalid='ignore'):
        if ((lat < -pi / 2) | (lat > pi / 2)).any():
            raise ValueError('-90 <= lat <= 90')

        if ((lon < -pi) | (lon > 2 * pi)).any():
            raise ValueError('-180 <= lat <= 360')

    if deg:
        return degrees(lat), degrees(lon), alt
    else:
        return lat, lon, alt  # radians


def ecef2enuv(u: float, v: float, w: float,
              lat0: float, lon0: float, deg: bool=True) -> Tuple[float, float, float]:
    """
    for VECTOR i.e. between two points

    input
    -----
    x,y,z  [meters] target ECEF location                             [0,Infinity)

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
             ell: Ellipsoid=None, deg: bool=True) -> Tuple[float, float, float]:
    """

    input
    -----
    x,y,z  [meters] target ECEF location                             [0,Infinity)
    Observer: lat0, lon0, h0 (altitude, meters)
    ell    reference ellipsoid
    deg    degrees input/output  (False: radians in/out)

    output:
    -------
    e,n,u   East, North, Up [m]

    """
    x0, y0, z0 = geodetic2ecef(lat0, lon0, h0, ell, deg=deg)

    return uvw2enu(x - x0, y - y0, z - z0, lat0, lon0, deg=deg)


def enu2uvw(east: float, north: float, up: float,
            lat0: float, lon0: float, deg: bool=True) -> Tuple[float, float, float]:
    if deg:
        lat0 = radians(lat0)
        lon0 = radians(lon0)

    t = cos(lat0) * up - sin(lat0) * north
    w = sin(lat0) * up + cos(lat0) * north

    u = cos(lon0) * t - sin(lon0) * east
    v = sin(lon0) * t + cos(lon0) * east

    return u, v, w


def uvw2enu(u: float, v: float, w: float,
            lat0: float, lon0: float, deg: bool=True) -> Tuple[float, float, float]:
    if deg:
        lat0 = radians(lat0)
        lon0 = radians(lon0)

    t = cos(lon0) * u + sin(lon0) * v
    East = -sin(lon0) * u + cos(lon0) * v
    Up = cos(lat0) * t + sin(lat0) * w
    North = -sin(lat0) * t + cos(lat0) * w

    return East, North, Up


def eci2geodetic(eci: np.ndarray, t: datetime) -> Tuple[float, float, float]:
    """
    convert ECI to geodetic coordinates

    inputs:

    eci/ecef: Nx3 vector of x,y,z triplets in the eci or ecef system [meters]
    t : length N vector of datetime OR greenwich sidereal time angle [radians].


    output
    ------
    lat,lon   (degrees/radians)
    alt  (meters)

    Note: Conversion is idealized: doesn't consider nutations, perterbations,
    etc. like the IAU-76/FK5 or IAU-2000/2006 model-based conversions
    from ECI to ECEF

    eci2geodetic() a.k.a. eci2lla()
    """
    ecef = np.atleast_2d(eci2ecef(eci, t))

    return np.asarray(ecef2geodetic(ecef[:, 0], ecef[:, 1], ecef[:, 2])).squeeze()


def enu2ecef(e1: float, n1: float, u1: float,
             lat0: float, lon0: float, h0: float,
             ell: Ellipsoid=None, deg: bool=True) -> Tuple[float, float, float]:
    """
    ENU to ECEF

    inputs:
     e1, n1, u1 (meters)   east, north, up
     observer: lat0, lon0, h0 (degrees/radians,degrees/radians, meters)
    ell    reference ellipsoid
    deg    degrees input/output  (False: radians in/out)


    output
    ------
    x,y,z  [meters] target ECEF location                         [0,Infinity)
    """
    x0, y0, z0 = geodetic2ecef(lat0, lon0, h0, ell, deg=deg)
    dx, dy, dz = enu2uvw(e1, n1, u1, lat0, lon0, deg=deg)

    return x0 + dx, y0 + dy, z0 + dz
