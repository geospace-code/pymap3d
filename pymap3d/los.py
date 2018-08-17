from typing import Tuple
import numpy as np
from math import pi
from .aer import aer2enu
from .ecef import enu2uvw, geodetic2ecef, Ellipsoid, ecef2geodetic


def lookAtSpheroid(lat0: float, lon0: float, h0: float, az: float, tilt: float,
                   ell=Ellipsoid(), deg: bool=True) -> Tuple[float, float, float]:
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
    if (np.asarray(h0) < 0).any():
        raise ValueError('Altitude \in  [0, Infinity)')

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
