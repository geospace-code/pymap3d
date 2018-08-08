from numpy import hypot, arctan2, degrees, radians, sin, cos
from typing import Tuple


def enu2aer(e: float, n: float, u: float, deg: bool=True) -> Tuple[float, float, float]:
    """
    Observer => Point

    input
    -----
    e,n,u  [meters]  East, north, up                                [0,Infinity)
    deg    degrees input/output  (False: radians in/out)

    output: AER
    ------
    azimuth, elevation (degrees/radians)                             [0,360),[0,90]
    slant range [meters]                                             [0,Infinity)
    """
    r = hypot(e, n)
    slantRange = hypot(r, u)
    elev = arctan2(u, r)
    az = arctan2(e, n) % (2 * arctan2(0, -1))
    if deg:
        return degrees(az), degrees(elev), slantRange
    else:
        return az, elev, slantRange  # radians


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


def aer2enu(az: float, el: float, srange: float, deg: bool=True) -> Tuple[float, float, float]:
    """
    azimuth, elevation (degrees/radians)                             [0,360),[0,90]
    slant range [meters]                                             [0,Infinity)
    deg    degrees input/output  (False: radians in/out)

    output:
    -------
    e,n,u   East, North, Up [m]

    """
    if deg:
        el = radians(el)
        az = radians(az)

    r = srange * cos(el)

    return r * sin(az), r * cos(az), srange * sin(el)


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
