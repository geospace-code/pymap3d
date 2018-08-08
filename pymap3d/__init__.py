#!/usr/bin/env python
# Copyright (c) 2014-2018 Michael Hirsch, Ph.D.

"""
Input/output: default units are METERS and DEGREES.
boolean deg=True means degrees

Most functions accept Numpy arrays of any shape

see tests/Test.py for example uses.
"""
from datetime import datetime
from typing import Tuple, Sequence, Union
import numpy as np
from math import pi
from .timeconv import str2dt  # noqa: F401
from .azelradec import radec2azel, azel2radec  # noqa: F401
from .aer import enu2aer, aer2enu, ned2aer, aer2ned  # noqa: F401
from .eci import eci2ecef, ecef2eci
from .datetime2hourangle import datetime2sidereal  # noqa: F401
from .ecef import geodetic2ecef, ecef2geodetic, EarthEllipsoid, ecef2enuv, enu2uvw, uvw2enu, ecef2enu, ecef2ned  # noqa: F401


def lookAtSpheroid(lat0: float, lon0: float, h0: float, az: float, tilt: float,
                   ell=EarthEllipsoid(), deg: bool=True) -> Tuple[float, float, float]:
    """
    Calculates line-of-sight intersection with Earth (or other ellipsoid) surface from above surface / orbit

    Args:
        lat0, lon0: latitude and longitude of starting point
        h0: altitude of starting point in meters
        az: azimuth angle of line-of-sight, clockwise from North
        tilt: tilt angle of line-of-sight with respect to local vertical (nadir = 0)
    Returns:
        lat, lon: latitude and longitude where the line-of-sight intersects with the Earth ellipsoid
        d: slant range in meters from the starting point to the intersect point
        Values will be NaN if the line of sight does not intersect.
    Algorithm based on https://medium.com/@stephenhartzell/satellite-line-of-sight-intersection-with-earth-d786b4a6a9b6 Stephen Hartzell
    """
    tilt = np.asarray(tilt)

    a = ell.a
    b = ell.a
    c = ell.b

    el = tilt - 90. if deg else tilt - pi / 2

    e, n, u = aer2enu(az, el, srange=1., deg=deg)  # fixed 1 km slant range
    u, v, w = enu2uvw(e, n, u, lat0, lon0, deg=deg)
    x, y, z = geodetic2ecef(lat0, lon0, h0, deg=deg)

    value = -a**2 * b**2 * w * z - a**2 * c**2 * v * y - b**2 * c**2 * u * x
    radical = (a**2 * b**2 * w**2 + a**2 * c**2 * v**2 - a**2 * v**2 * z**2 + 2 * a**2 * v * w * y * z -
               a**2 * w**2 * y**2 + b**2 * c**2 * u**2 - b**2 * u**2 * z**2 + 2 * b**2 * u * w * x * z -
               b**2 * w**2 * x**2 - c**2 * u**2 * y**2 + 2 * c**2 * u * v * x * y - c**2 * v**2 * x**2)

    magnitude = a**2 * b**2 * w**2 + a**2 * c**2 * v**2 + b**2 * c**2 * u**2

#   Return nan if radical < 0 or d < 0 because LOS vector does not point towards Earth
    with np.errstate(invalid='ignore'):
        d = np.where(radical > 0,
                     (value - a * b * c * np.sqrt(radical)) / magnitude,
                     np.nan)
        d[d < 0] = np.nan

    lat, lon, _ = ecef2geodetic(x + d * u, y + d * v, z + d * w, deg=deg)

    return lat, lon, d


# alias
los_intersect = lookAtSpheroid
# %% to AER (azimuth, elevation, range)


def ecef2aer(x: float, y: float, z: float,
             lat0: float, lon0: float, h0: float,
             ell=None, deg: bool=True) -> Tuple[float, float, float]:
    """
    Observer => Point

    input:
    -----
    x,y,z  [meters] target ECEF location                             [0,Infinity)
    lat0, lon0 (degrees/radians)  Observer coordinates on ellipsoid  [-90,90],[-180,180]
    h0     [meters]                observer altitude                 [0,Infinity)
    ell    reference ellipsoid
    deg    degrees input/output  (False: radians in/out)

    output: AER
    ------
    azimuth, elevation (degrees/radians)                             [0,360),[0,90]
    slant range [meters]                                             [0,Infinity)
    """
    xEast, yNorth, zUp = ecef2enu(x, y, z, lat0, lon0, h0, ell, deg=deg)

    return enu2aer(xEast, yNorth, zUp, deg=deg)


def geodetic2aer(lat: float, lon: float, h: float,
                 lat0: float, lon0: float, h0: float,
                 ell=None, deg: bool=True) -> Tuple[float, float, float]:
    """
    Observer => Point

    input:
    -----
    Target:   lat, lon, h (altitude, meters)
    Observer: lat0, lon0, h0 (altitude, meters)
    ell    reference ellipsoid
    deg    degrees input/output  (False: radians in/out)

    output: AER
    ------
    azimuth, elevation (degrees/radians)
    slant range [meters]
    """
    e, n, u = geodetic2enu(lat, lon, h, lat0, lon0, h0, ell, deg=deg)

    return enu2aer(e, n, u, deg=deg)

# %% to ECEF


def aer2ecef(az: float, el: float, srange: float,
             lat0: float, lon0: float, alt0: float,
             ell=None, deg: bool=True) -> Tuple[float, float, float]:
    """
    convert target azimuth, elevation, range (meters) from observer at lat0,lon0,alt0 to ECEF coordinates.

     Input:
     -----
    azimuth, elevation (degrees/radians)                             [0,360),[0,90]
    slant range [meters]                                             [0,Infinity)
    Observer: lat0, lon0, h0 (altitude, meters)
    ell    reference ellipsoid
    deg    degrees input/output  (False: radians in/out)

     output: ECEF x,y,z  [meters]

    if you specify NaN for srange, return value z will be NaN
    """
    # Origin of the local system in geocentric coordinates.
    x0, y0, z0 = geodetic2ecef(lat0, lon0, alt0, ell, deg=deg)
    # Convert Local Spherical AER to ENU
    e1, n1, u1 = aer2enu(az, el, srange, deg=deg)
    # Rotating ENU to ECEF
    dx, dy, dz = enu2uvw(e1, n1, u1, lat0, lon0, deg=deg)
    # Origin + offset from origin equals position in ECEF
    return x0 + dx, y0 + dy, z0 + dz


def enu2ecef(e1: float, n1: float, u1: float,
             lat0: float, lon0: float, h0: float,
             ell=None, deg: bool=True) -> Tuple[float, float, float]:
    """
    Observer => Point

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


def ned2ecef(n: float, e: float, d: float,
             lat0: float, lon0: float, h0: float,
             ell=None, deg: bool=True) -> Tuple[float, float, float]:
    """
    Observer => Point

    input
    -----
    n,e,d  [meters]  North,east, down                                [0,Infinity)
    Observer: lat0, lon0, h0 (altitude, meters)
    ell    reference ellipsoid
    deg    degrees input/output  (False: radians in/out)

    output:
    ------
    ECEF  x,y,z  (meters)
    """
    return enu2ecef(e, n, -d, lat0, lon0, h0, ell, deg=deg)

# %% ECI


def aer2eci(az: float, el: float, srange: float,
            lat0: float, lon0: float, h0: float, t: datetime,
            ell=None, deg: bool=True) -> np.ndarray:
    """

    input
    -----
    azimuth, elevation (degrees/radians)                             [0,360),[0,90]
    slant range [meters]                                             [0,Infinity)
    Observer: lat0, lon0, h0 (altitude, meters)
    ell    reference ellipsoid
    deg    degrees input/output  (False: radians in/out)
    t  datetime.datetime of obseration

    output
    ------
    eci  x,y,z (meters)
    """
    x, y, z = aer2ecef(az, el, srange, lat0, lon0, h0, ell, deg)

    return ecef2eci(np.column_stack((x, y, z)), t)


def eci2aer(eci: Sequence[float], lat0: float, lon0: float, h0: float,
            t: Union[datetime, Sequence[datetime]]) -> Tuple[float, float, float]:
    """
    Observer => Point

    input
    -----
    eci [meters] Nx3 target ECI location (x,y,z)                    [0,Infinity)
    lat0, lon0 (degrees/radians)  Observer coordinates on ellipsoid [-90,90],[-180,180]
    h0     [meters]                observer altitude                [0,Infinity)
    t  time (datetime.datetime)   time of obsevation (UTC)


    output: AER
    ------
    azimuth, elevation (degrees/radians)                             [0,360),[0,90]
    slant range [meters]                                             [0,Infinity)
    """
    ecef = eci2ecef(eci, t)

    return ecef2aer(ecef[:, 0], ecef[:, 1], ecef[:, 2], lat0, lon0, h0)

# %% to ENU


def geodetic2enu(lat: float, lon: float, h: float,
                 lat0: float, lon0: float, h0: float,
                 ell=None, deg: bool=True) -> Tuple[float, float, float]:
    """
    input
    -----
    target: lat,lon, h
    Observer: lat0, lon0, h0 (altitude, meters)
    ell    reference ellipsoid
    deg    degrees input/output  (False: radians in/out)


    output:
    -------
    e,n,u   East, North, Up [m]
    """
    x1, y1, z1 = geodetic2ecef(lat, lon, h, ell, deg=deg)
    x2, y2, z2 = geodetic2ecef(lat0, lon0, h0, ell, deg=deg)
    dx = x1 - x2
    dy = y1 - y2
    dz = z1 - z2

    return uvw2enu(dx, dy, dz, lat0, lon0, deg=deg)
# %% to geodetic


def aer2geodetic(az: float, el: float, srange: float,
                 lat0: float, lon0: float, h0: float, deg: bool=True) -> Tuple[float, float, float]:
    """
    Input:
    -----
    az,el (degrees/radians)
    srange[meters]

    Observer: lat0,lon0 [degrees]
              altitude h0 [meters]

    deg :   degrees input/output  (False: radians in/out)

    output:
    WGS84 lat,lon [degrees]  h0altitude above spheroid  [meters]

    """
    x, y, z = aer2ecef(az, el, srange, lat0, lon0, h0, deg=deg)

    return ecef2geodetic(x, y, z, deg=deg)


def eci2geodetic(eci: Sequence[float], t: datetime) -> Tuple[float, float, float]:
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
    """

    """ a.k.a. eci2lla() """
    ecef = eci2ecef(eci, t)

    return ecef2geodetic(ecef[:, 0], ecef[:, 1], ecef[:, 2])


def enu2geodetic(e: float, n: float, u: float,
                 lat0: float, lon0: float, h0: float,
                 ell=None, deg: bool=True) -> Tuple[float, float, float]:
    """

    input
    -----
    e,n,u   East, North, Up [m]
    Observer: lat0, lon0, h0 (altitude, meters)
    ell    reference ellipsoid
    deg    degrees input/output  (False: radians in/out)


    output:
    -------
    target: lat,lon, h  (degrees/radians,degrees/radians, meters)

    """

    x, y, z = enu2ecef(e, n, u, lat0, lon0, h0, ell, deg=deg)

    return ecef2geodetic(x, y, z, ell, deg=deg)


def ned2geodetic(n: float, e: float, d: float,
                 lat0: float, lon0: float, h0: float,
                 ell=None, deg: bool=True) -> Tuple[float, float, float]:
    """

    input
    -----
    n,e,d   North, east, down (meters)
    Observer: lat0, lon0, h0 (altitude, meters)
    ell    reference ellipsoid
    deg    degrees input/output  (False: radians in/out)


    output:
    -------
    target: lat,lon, h  (degrees/radians,degrees/radians, meters)

    """
    x, y, z = enu2ecef(e, n, -d, lat0, lon0, h0, ell, deg=deg)

    return ecef2geodetic(x, y, z, ell, deg=deg)
# %% to NED


def ecef2nedv(u, v, w, lat0, lon0, deg=True) -> Tuple[float, float, float]:
    """
    for VECTOR between two points
    """
    e, n, u = ecef2enuv(u, v, w, lat0, lon0, deg=deg)

    return n, e, -u


def geodetic2ned(lat, lon, h, lat0, lon0, h0, ell=None, deg=True) -> Tuple[float, float, float]:
    """
    input
    -----
    target: lat,lon  (degrees/radians)
    h (altitude, meters)
    Observer: lat0, lon0 (degrees/radians)
    h0 (altitude, meters)
    ell    reference ellipsoid
    deg    degrees input/output  (False: radians in/out)


    output:
    -------
    n,e,d  North,east, down [m]
    """
    e, n, u = geodetic2enu(lat, lon, h, lat0, lon0, h0, ell, deg=deg)

    return n, e, -u
