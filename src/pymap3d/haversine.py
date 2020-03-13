"""
Compute angular separation in the sky using haversine

Note:
decimal points on constants made 0 difference in `%timeit` execution time

The Meeus algorithm is about 9.5% faster than Astropy/Vicenty on my PC,
and gives virtually identical result
within double precision arithmetic limitations
"""

import typing

try:
    from numpy import cos, arcsin, sqrt, radians, degrees
except ImportError:
    from math import cos, sqrt, radians, degrees, asin as arcsin
try:
    from astropy.coordinates.angle_utilities import angular_separation
except ImportError:
    angular_separation = None

__all__ = ["anglesep", "anglesep_meeus", "haversine"]

if typing.TYPE_CHECKING:
    from numpy import ndarray


def anglesep_meeus(lon0: "ndarray", lat0: "ndarray", lon1: "ndarray", lat1: "ndarray", deg: bool = True) -> "ndarray":
    """
    Parameters
    ----------

    lon0 : "ndarray"
        longitude of first point
    lat0 : "ndarray"
        latitude of first point
    lon1 : "ndarray"
        longitude of second point
    lat1 : "ndarray"
        latitude of second point
    deg : bool, optional
          degrees input/output  (False: radians in/out)

    Returns
    -------

    sep_rad : "ndarray"
        angular separation


    Meeus p. 109

    from "Astronomical Algorithms" by Jean Meeus Ch. 16 p. 111 (16.5)
    gives angular distance in degrees between two rightAscension,Declination
    points in the sky.  Neglecting atmospheric effects, of course.

    Meeus haversine method is stable all the way to exactly 0 deg.

    either the arrays must be the same size, or one of them must be a scalar
    """

    if deg:
        lon0 = radians(lon0)
        lat0 = radians(lat0)
        lon1 = radians(lon1)
        lat1 = radians(lat1)

    sep_rad = 2 * arcsin(sqrt(haversine(lat0 - lat1) + cos(lat0) * cos(lat1) * haversine(lon0 - lon1)))

    return degrees(sep_rad) if deg else sep_rad


def anglesep(lon0: "ndarray", lat0: "ndarray", lon1: "ndarray", lat1: "ndarray", deg: bool = True) -> "ndarray":
    """
    Parameters
    ----------

    lon0 : "ndarray"
        longitude of first point
    lat0 : "ndarray"
        latitude of first point
    lon1 : "ndarray"
        longitude of second point
    lat1 : "ndarray"
        latitude of second point
    deg : bool, optional
          degrees input/output  (False: radians in/out)

    Returns
    -------

    sep_rad : "ndarray"
        angular separation

    For reference, this is from astropy astropy/coordinates/angle_utilities.py
    Angular separation between two points on a sphere.
    """
    if angular_separation is None:
        raise ImportError("angledist requires AstroPy. Try angledis_meeus")

    if deg:
        lon0 = radians(lon0)
        lat0 = radians(lat0)
        lon1 = radians(lon1)
        lat1 = radians(lat1)

    sep_rad = angular_separation(lon0, lat0, lon1, lat1)

    return degrees(sep_rad) if deg else sep_rad


def haversine(theta: "ndarray") -> "ndarray":
    """
    Compute haversine

    Parameters
    ----------

    theta : "ndarray"
        angle (radians)

    Results
    -------

    htheta : "ndarray"
        haversine of `theta`

    https://en.wikipedia.org/wiki/Haversine
    Meeus p. 111
    """
    return (1 - cos(theta)) / 2.0
