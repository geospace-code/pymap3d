from .enu import geodetic2enu, aer2enu, enu2aer
from .ecef import ecef2geodetic, ecef2enuv, ecef2enu, enu2ecef
from typing import Tuple


def aer2ned(az: float, elev: float, slantRange: float, deg: bool=True) -> Tuple[float, float, float]:
    """
    input
    -----
    azimuth, elevation (degrees/radians)                             [0,360),[0,90]
    slant range [meters]                                             [0,Infinity)
    deg    degrees input/output  (False: radians in/out)

    output:
    -------
    n,e,d  North,east, down [m]
    """
    e, n, u = aer2enu(az, elev, slantRange, deg=deg)

    return n, e, -u


def ned2aer(n: float, e: float, d: float, deg: bool=True) -> Tuple[float, float, float]:
    """
    Observer => Point

    input
    -----
    n,e,d  [meters]  North,east, down                                [0,Infinity)
    deg    degrees input/output  (False: radians in/out)

    output: AER
    ------
    azimuth, elevation (degrees/radians)                             [0,360),[0,90]
    slant range [meters]                                             [0,Infinity)
    """
    return enu2aer(e, n, -d, deg=deg)


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


def ecef2nedv(u, v, w, lat0, lon0, deg=True) -> Tuple[float, float, float]:
    """
    for VECTOR between two points
    """
    e, n, u = ecef2enuv(u, v, w, lat0, lon0, deg=deg)

    return n, e, -u
