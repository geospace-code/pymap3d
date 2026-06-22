"""compute radii of curvature for an ellipsoid"""

from __future__ import annotations

from ._typing import FloatLike, FloatArray

from .ellipsoid import Ellipsoid, resolve_ellipsoid
from .mathfun import cos, sin, sqrt, radians

__all__ = ["parallel", "meridian", "transverse", "geocentric_radius"]


def geocentric_radius(geodetic_lat: FloatLike, ell: Ellipsoid | None = None, deg: bool = True):
    """
    compute geocentric radius at geodetic latitude

    https://en.wikipedia.org/wiki/Earth_radius#Geocentric_radius
    """

    ell = resolve_ellipsoid(ell)

    if deg:
        geodetic_lat = radians(geodetic_lat)

    return sqrt(
        (
            (ell.semimajor_axis**2 * cos(geodetic_lat)) ** 2
            + (ell.semiminor_axis**2 * sin(geodetic_lat)) ** 2
        )
        / (
            (ell.semimajor_axis * cos(geodetic_lat)) ** 2
            + (ell.semiminor_axis * sin(geodetic_lat)) ** 2
        )
    )


def parallel(lat: FloatLike | FloatArray, ell: Ellipsoid | None = None, deg: bool = True):
    """
    computes the radius of the small circle encompassing the globe at the specified latitude

    like Matlab rcurve('parallel', ...)

    Parameters
    ----------
    lat : float
        geodetic latitude (degrees)
    ell : Ellipsoid, optional
          reference ellipsoid
    deg : bool, optional
          degrees input/output  (False: radians in/out)

    Returns
    -------
    radius: float
        radius of ellipsoid (meters)
    """

    if deg:
        lat = radians(lat)

    return cos(lat) * transverse(lat, ell, deg=False)


def meridian(lat, ell: Ellipsoid | None = None, deg: bool = True):
    """computes the meridional radius of curvature for the ellipsoid

    like Matlab rcurve('meridian', ...)

    Parameters
    ----------
    lat : array-like float
        geodetic latitude (degrees)
    ell : Ellipsoid, optional
          reference ellipsoid
    deg : bool, optional
          degrees input/output  (False: radians in/out)

    Returns
    -------
    radius: array-like float
        radius of ellipsoid
    """

    ell = resolve_ellipsoid(ell)

    if deg:
        lat = radians(lat)

    f1 = ell.semimajor_axis * (1 - ell.eccentricity**2)
    f2 = 1 - (ell.eccentricity * sin(lat)) ** 2
    return f1 / sqrt(f2**3)


def transverse(lat: FloatLike | FloatArray, ell: Ellipsoid | None = None, deg: bool = True):
    """computes the radius of the curve formed by a plane
    intersecting the ellipsoid at the latitude which is
    normal to the surface of the ellipsoid

    like Matlab rcurve('transverse', ...)

    Parameters
    ----------
    lat : array-like float
        latitude (degrees)
    ell : Ellipsoid, optional
          reference ellipsoid
    deg : bool, optional
          degrees input/output  (False: radians in/out)

    Returns
    -------
    radius: array-like float
        radius of ellipsoid (meters)
    """

    ell = resolve_ellipsoid(ell)

    if deg:
        lat = radians(lat)

    return ell.semimajor_axis / sqrt(1 - (ell.eccentricity * sin(lat)) ** 2)
