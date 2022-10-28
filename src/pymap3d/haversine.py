"""
Compute angular separation in the sky using haversine

Note:
decimal points on constants made 0 difference in `%timeit` execution time

The Meeus algorithm is about 9.5% faster than Astropy/Vicenty on my PC,
and gives virtually identical result
within double precision arithmetic limitations
"""
from __future__ import annotations

from typing import Any, Sequence, overload

try:
    from astropy.coordinates.angle_utilities import angular_separation
except ImportError:
    pass

try:
    from numpy import asarray
    from numpy.typing import NDArray
except ImportError:
    pass

from ._types import ArrayLike
from .mathfun import asin, cos, degrees, radians, sqrt

__all__ = ["anglesep", "anglesep_meeus", "haversine"]


@overload
def anglesep_meeus(
    lon0: float,
    lat0: float,
    lon1: float,
    lat1: float,
    deg: bool = True,
) -> float:
    pass


@overload
def anglesep_meeus(
    lon0: ArrayLike,
    lat0: ArrayLike,
    lon1: ArrayLike,
    lat1: ArrayLike,
    deg: bool = True,
) -> NDArray[Any]:
    pass


def anglesep_meeus(
    lon0: float | ArrayLike,
    lat0: float | ArrayLike,
    lon1: float | ArrayLike,
    lat1: float | ArrayLike,
    deg: bool = True,
) -> float | NDArray[Any]:
    """
    Parameters
    ----------

    lon0 : float
        longitude of first point
    lat0 : float
        latitude of first point
    lon1 : float
        longitude of second point
    lat1 : float
        latitude of second point
    deg : bool, optional
          degrees input/output  (False: radians in/out)

    Returns
    -------

    sep_rad : float
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
    else:
        try:
            lon0 = asarray(lon0)
            lat0 = asarray(lat0)
            lon1 = asarray(lon1)
            lat1 = asarray(lat1)
        except NameError:
            pass
    assert (
        not isinstance(lon0, Sequence)
        and not isinstance(lat0, Sequence)
        and not isinstance(lon1, Sequence)
        and not isinstance(lat1, Sequence)
    )

    sep_rad = 2 * asin(
        sqrt(haversine(lat0 - lat1) + cos(lat0) * cos(lat1) * haversine(lon0 - lon1))
    )

    return degrees(sep_rad) if deg else sep_rad  # type: ignore[no-any-return]


@overload
def anglesep(
    lon0: float,
    lat0: float,
    lon1: float,
    lat1: float,
    deg: bool = True,
) -> float:
    pass


@overload
def anglesep(
    lon0: ArrayLike,
    lat0: ArrayLike,
    lon1: ArrayLike,
    lat1: ArrayLike,
    deg: bool = True,
) -> NDArray[Any]:
    pass


def anglesep(
    lon0: float | ArrayLike,
    lat0: float | ArrayLike,
    lon1: float | ArrayLike,
    lat1: float | ArrayLike,
    deg: bool = True,
) -> float | NDArray[Any]:
    """
    Parameters
    ----------

    lon0 : float
        longitude of first point
    lat0 : float
        latitude of first point
    lon1 : float
        longitude of second point
    lat1 : float
        latitude of second point
    deg : bool, optional
          degrees input/output  (False: radians in/out)

    Returns
    -------

    sep_rad : float
        angular separation

    For reference, this is from astropy astropy/coordinates/angle_utilities.py
    Angular separation between two points on a sphere.
    """

    if deg:
        lon0 = radians(lon0)
        lat0 = radians(lat0)
        lon1 = radians(lon1)
        lat1 = radians(lat1)
    else:
        try:
            lon0 = asarray(lon0)
            lat0 = asarray(lat0)
            lon1 = asarray(lon1)
            lat1 = asarray(lat1)
        except NameError:
            pass
    assert (
        not isinstance(lon0, Sequence)
        and not isinstance(lat0, Sequence)
        and not isinstance(lon1, Sequence)
        and not isinstance(lat1, Sequence)
    )

    try:
        sep_rad = angular_separation(lon0, lat0, lon1, lat1)
    except NameError:
        sep_rad = anglesep_meeus(lon0, lat0, lon1, lat1, deg=False)  # type: ignore[arg-type]

    return degrees(sep_rad) if deg else sep_rad  # type: ignore[no-any-return]


@overload
def haversine(theta: float) -> float:
    pass


@overload
def haversine(theta: ArrayLike) -> NDArray[Any]:
    pass


def haversine(theta: float | ArrayLike) -> float | NDArray[Any]:
    """
    Compute haversine

    Parameters
    ----------

    theta : float
        angle (radians)

    Results
    -------

    htheta : float
        haversine of `theta`

    https://en.wikipedia.org/wiki/Haversine
    Meeus p. 111
    """
    return (1 - cos(theta)) / 2.0
