"""compute radii of curvature for an ellipsoid"""

import typing

try:
    from numpy import radians, sin, cos, sqrt
except ImportError:
    from math import radians, sin, cos, sqrt
from .ellipsoid import Ellipsoid

__all__ = ["rcurve_parallel", "rcurve_meridian", "rcurve_transverse"]

if typing.TYPE_CHECKING:
    from numpy import ndarray


def rcurve_parallel(lat: "ndarray", ell: Ellipsoid = None, deg: bool = True) -> "ndarray":
    """
    computes the radius of the small circle encompassing the globe at the specified latitude

    Parameters
    ----------
    lat : "ndarray"
        geodetic latitude (degrees)
    ell : Ellipsoid, optional
          reference ellipsoid
    deg : bool, optional
          degrees input/output  (False: radians in/out)

    Returns
    -------
    radius: "ndarray"
        radius of ellipsoid
    """

    if deg:
        lat = radians(lat)

    return cos(lat) * rcurve_transverse(lat, ell, deg=False)


def rcurve_meridian(lat: "ndarray", ell: Ellipsoid = None, deg: bool = True) -> "ndarray":
    """computes the meridional radius of curvature for the ellipsoid

    Parameters
    ----------
    lat : "ndarray"
        geodetic latitude (degrees)
    ell : Ellipsoid, optional
          reference ellipsoid
    deg : bool, optional
          degrees input/output  (False: radians in/out)

    Returns
    -------
    radius: "ndarray"
        radius of ellipsoid
    """

    if ell is None:
        ell = Ellipsoid()
    if deg:
        lat = radians(lat)

    f1 = ell.semimajor_axis * (1 - ell.eccentricity ** 2)
    f2 = 1 - (ell.eccentricity * sin(lat)) ** 2
    return f1 / sqrt(f2 ** 3)


def rcurve_transverse(lat: "ndarray", ell: Ellipsoid = None, deg: bool = True) -> "ndarray":
    """computes the radius of the curve formed by a plane
    intersecting the ellipsoid at the latitude which is
    normal to the surface of the ellipsoid

    Parameters
    ----------
    lat : "ndarray"
        geodetic latitude (degrees)
    ell : Ellipsoid, optional
          reference ellipsoid
    deg : bool, optional
          degrees input/output  (False: radians in/out)

    Returns
    -------
    radius: "ndarray"
        radius of ellipsoid
    """

    if ell is None:
        ell = Ellipsoid()
    if deg:
        lat = radians(lat)

    return ell.semimajor_axis / sqrt(1 - (ell.eccentricity * sin(lat)) ** 2)
