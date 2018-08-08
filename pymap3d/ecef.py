from numpy import radians, sin, cos, tan, allclose, hypot, degrees, arctan2, sqrt
from copy import deepcopy
from typing import Tuple


class EarthEllipsoid:
    """generate reference ellipsoid"""

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


def get_radius_normal(lat_radians: float, ell=None) -> float:
    """ Compute normal radius of planetary body"""
    if ell is None:
        ell = EarthEllipsoid()

    a = ell.a
    b = ell.b

    return a**2 / sqrt(a**2 * cos(lat_radians)**2 + b**2 * sin(lat_radians)**2)


def geodetic2ecef(lat: float, lon: float, alt: float,
                  ell=None, deg: bool=True) -> Tuple[float, float, float]:
    """
    Observer => Point

    input:
    -----
    Target:   lat, lon, h (altitude, meters)
    Observer: lat0, lon0, h0 (altitude, meters)
    ell    reference ellipsoid
    deg    degrees input/output  (False: radians in/out)

    output: ECEF x,y,z (meters)
    """
    if ell is None:
        ell = EarthEllipsoid()

    if deg:
        lat = radians(lat)
        lon = radians(lon)
    # radius of curvature of the prime vertical section
    N = get_radius_normal(lat, ell)
    # Compute cartesian (geocentric) coordinates given  (curvilinear) geodetic
    # coordinates.
    x = (N + alt) * cos(lat) * cos(lon)
    y = (N + alt) * cos(lat) * sin(lon)
    z = (N * (ell.b / ell.a)**2 + alt) * sin(lat)

    return x, y, z


def ecef2geodetic(x: float, y: float, z: float,
                  ell=None, deg: bool=True) -> Tuple[float, float, float]:
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
        ell = EarthEllipsoid()

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
             ell=None, deg: bool=True) -> Tuple[float, float, float]:
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


def ecef2ned(x, y, z, lat0, lon0, h0, ell=None, deg=True) -> Tuple[float, float, float]:
    """
    input
    -----
    x,y,z  [meters] target ECEF location                             [0,Infinity)
    Observer: lat0, lon0, h0 (altitude, meters)
    ell    reference ellipsoid
    deg    degrees input/output  (False: radians in/out)

    output:
    -------
    n,e,d  North,east, down [m]

    """
    e, n, u = ecef2enu(x, y, z, lat0, lon0, h0, ell, deg=deg)

    return n, e, -u


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
