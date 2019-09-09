""" Line of sight intersection of space observer to ellipsoid """
from typing import Tuple
from math import pi, nan, sqrt

from .aer import aer2enu
from .ecef import enu2uvw, geodetic2ecef, ecef2geodetic
from .ellipsoid import Ellipsoid
try:
    import numpy
except ImportError:
    numpy = None

__all__ = ['lookAtSpheroid']


def lookAtSpheroid(lat0: float, lon0: float, h0: float, az: float, tilt: float,
                   ell: Ellipsoid = None, deg: bool = True) -> Tuple[float, float, float]:
    if numpy is not None:
        fun = numpy.vectorize(lookAtSpheroid_point)
        return fun(lat0, lon0, h0, az, tilt, ell, deg)
    else:
        return lookAtSpheroid_point(lat0, lon0, h0, az, tilt, ell, deg)


def lookAtSpheroid_point(lat0: float, lon0: float, h0: float, az: float, tilt: float,
                         ell: Ellipsoid = None, deg: bool = True) -> Tuple[float, float, float]:
    """
    Calculates line-of-sight intersection with Earth (or other ellipsoid) surface from above surface / orbit

    Parameters
    ----------

    lat0 : float
           observer geodetic latitude
    lon0 : float
           observer geodetic longitude
    h0 : float
        observer altitude (meters)  Must be non-negative since this function doesn't consider terrain
    az : float or numpy.ndarray of float
        azimuth angle of line-of-sight, clockwise from North
    tilt : float or numpy.ndarray of float
        tilt angle of line-of-sight with respect to local vertical (nadir = 0)
    ell : Ellipsoid, optional
          reference ellipsoid
    deg : bool, optional
          degrees input/output  (False: radians in/out)

    Results
    -------

    lat0 : float or numpy.ndarray of float
           geodetic latitude where the line-of-sight intersects with the Earth ellipsoid
    lon0 : float or numpy.ndarray of float
           geodetic longitude where the line-of-sight intersects with the Earth ellipsoid
    d : float or numpy.ndarray of float
        slant range (meters) from starting point to intersect point

    Values will be NaN if the line of sight does not intersect.

    Algorithm based on https://medium.com/@stephenhartzell/satellite-line-of-sight-intersection-with-earth-d786b4a6a9b6 Stephen Hartzell
    """

    if h0 < 0:
        raise ValueError('Intersection calculation requires altitude  [0, Infinity)')

    if ell is None:
        ell = Ellipsoid()

    a = ell.semimajor_axis
    b = ell.semimajor_axis
    c = ell.semiminor_axis

    el = tilt - 90. if deg else tilt - pi / 2

    e, n, u = aer2enu(az, el, srange=1., deg=deg)  # fixed 1 km slant range
    u, v, w = enu2uvw(e, n, u, lat0, lon0, deg=deg)
    x, y, z = geodetic2ecef(lat0, lon0, h0, deg=deg)

    value = -a**2 * b**2 * w * z - a**2 * c**2 * v * y - b**2 * c**2 * u * x
    radical = (a**2 * b**2 * w**2 + a**2 * c**2 * v**2 - a**2 * v**2 * z**2 + 2 * a**2 * v * w * y * z -
               a**2 * w**2 * y**2 + b**2 * c**2 * u**2 - b**2 * u**2 * z**2 + 2 * b**2 * u * w * x * z -
               b**2 * w**2 * x**2 - c**2 * u**2 * y**2 + 2 * c**2 * u * v * x * y - c**2 * v**2 * x**2)

    magnitude = a**2 * b**2 * w**2 + a**2 * c**2 * v**2 + b**2 * c**2 * u**2

# %%   Return nan if radical < 0 or d < 0 because LOS vector does not point towards Earth
    if radical > 0:
        d = (value - a * b * c * sqrt(radical)) / magnitude
    else:
        d = nan

    if d < 0:
        d = nan
# %% cartesian to ellipsodal
    lat, lon, _ = ecef2geodetic(x + d * u, y + d * v, z + d * w, deg=deg)

    return lat, lon, d
