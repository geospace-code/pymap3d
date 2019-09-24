"""compute radii of curvature for an ellipsoid"""
try:
    from numpy import radians, sin, cos, sqrt
except ImportError:
    from math import radians, sin, cos, sqrt
from .ellipsoid import Ellipsoid

__all__ = ["rcurve_parallel", "rcurve_meridian", "rcurve_transverse"]


def rcurve_parallel(lat: float, ell: Ellipsoid = None, deg: bool = True) -> float:
    """
    computes the radius of the small circle encompassing the globe at the specified latitude
    """

    if deg:
        lat = radians(lat)

    return cos(lat) * rcurve_transverse(lat, ell, deg=False)


def rcurve_meridian(lat: float, ell: Ellipsoid = None, deg: bool = True) -> float:
    """computes the meridional radius of curvature for the ellipsoid"""

    if ell is None:
        ell = Ellipsoid()
    if deg:
        lat = radians(lat)

    f1 = ell.semimajor_axis * (1 - ell.eccentricity ** 2)
    f2 = 1 - (ell.eccentricity * sin(lat)) ** 2
    return f1 / sqrt(f2 ** 3)


def rcurve_transverse(lat: float, ell: Ellipsoid = None, deg: bool = True) -> float:
    """computes the radius of the curve formed by a plane
    intersecting the ellipsoid at the latitude which is
    normal to the surface of the ellipsoid
    """

    if ell is None:
        ell = Ellipsoid()
    if deg:
        lat = radians(lat)

    return ell.semimajor_axis / sqrt(1 - (ell.eccentricity * sin(lat)) ** 2)
