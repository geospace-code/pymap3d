""" compute radii of auxiliary spheres"""

import typing

try:
    from numpy import radians, sin, cos, log, sqrt, degrees, asarray
except ImportError:
    from math import radians, sin, cos, log, sqrt, degrees

    asarray = None

from .ellipsoid import Ellipsoid
from .rcurve import rcurve_meridian, rcurve_transverse
from .vincenty import vdist

try:
    from numpy.typing import ArrayLike
except ImportError:
    ArrayLike = typing.Any

__all__ = [
    "rsphere_eqavol",
    "rsphere_authalic",
    "rsphere_rectifying",
    "rsphere_euler",
    "rsphere_curve",
    "rsphere_triaxial",
    "rsphere_biaxial",
]


def rsphere_eqavol(ell: Ellipsoid = None) -> float:
    """computes the radius of the sphere with equal volume as the ellipsoid

    Parameters
    ----------
    ell : Ellipsoid, optional
          reference ellipsoid

    Returns
    -------
    radius: float
        radius of sphere
    """
    if ell is None:
        ell = Ellipsoid()

    f = ell.flattening

    return ell.semimajor_axis * (1 - f / 3 - f ** 2 / 9)


def rsphere_authalic(ell: Ellipsoid = None) -> float:
    """computes the radius of the sphere with equal surface area as the ellipsoid

    Parameters
    ----------
    ell : Ellipsoid, optional
          reference ellipsoid

    Returns
    -------
    radius: float
        radius of sphere
    """
    if ell is None:
        ell = Ellipsoid()

    e = ell.eccentricity

    if e > 0:
        f1 = ell.semimajor_axis ** 2 / 2
        f2 = (1 - e ** 2) / (2 * e)
        f3 = log((1 + e) / (1 - e))
        return sqrt(f1 * (1 + f2 * f3))
    else:
        return ell.semimajor_axis


def rsphere_rectifying(ell: Ellipsoid = None) -> float:
    """computes the radius of the sphere with equal meridional distances as the ellipsoid

    Parameters
    ----------
    ell : Ellipsoid, optional
          reference ellipsoid

    Returns
    -------
    radius: float
        radius of sphere
    """
    if ell is None:
        ell = Ellipsoid()
    return ((ell.semimajor_axis ** (3 / 2) + ell.semiminor_axis ** (3 / 2)) / 2) ** (2 / 3)


def rsphere_euler(
    lat1: ArrayLike, lon1: ArrayLike, lat2: ArrayLike, lon2: ArrayLike, ell: Ellipsoid = None, deg: bool = True
) -> ArrayLike:
    """computes the Euler radii of curvature at the midpoint of the
     great circle arc defined by the endpoints (lat1,lon1) and (lat2,lon2)

    Parameters
    ----------
    lat1, lat2 : ArrayLike
        geodetic latitudes (degrees)
    lon1, lon2 : ArrayLike
       geodetic longitudes (degrees)
    ell : Ellipsoid, optional
          reference ellipsoid
    deg : bool, optional
          degrees input/output  (False: radians in/out)

    Returns
    -------
    radius: ArrayLike
        radius of sphere
    """
    if not deg:
        lat1, lon1, lat2, lon2 = degrees(lat1), degrees(lon1), degrees(lat2), degrees(lon2)
    if asarray is not None:
        lat1, lat2 = asarray(lat1), asarray(lat2)

    latmid = lat1 + (lat2 - lat1) / 2  # compute the midpoint

    # compute azimuth
    az = vdist(lat1, lon1, lat2, lon2, ell=ell)[1]

    #   compute meridional and transverse radii of curvature
    rho = rcurve_meridian(latmid, ell, deg=True)
    nu = rcurve_transverse(latmid, ell, deg=True)

    az = radians(az)
    den = rho * sin(az) ** 2 + nu * cos(az) ** 2

    #  compute radius of the arc from point 1 to point 2
    return rho * nu / den


def rsphere_curve(lat: ArrayLike, ell: Ellipsoid = None, deg: bool = True, method: str = "mean") -> ArrayLike:
    """computes the arithmetic average of the transverse and meridional
    radii of curvature at a specified latitude point

    Parameters
    ----------
    lat1 : ArrayLike
        geodetic latitudes (degrees)
    ell : Ellipsoid, optional
          reference ellipsoid
    method: str, optional
        "mean" or "norm"
    deg : bool, optional
          degrees input/output  (False: radians in/out)

    Returns
    -------
    radius: ArrayLike
        radius of sphere
    """

    if deg:
        lat = radians(lat)

    rho = rcurve_meridian(lat, ell, deg=False)
    nu = rcurve_transverse(lat, ell, deg=False)

    if method == "mean":
        return (rho + nu) / 2
    elif method == "norm":
        return sqrt(rho * nu)
    else:
        raise Exception("pymap3d.rsphere.curve: method must be mean or norm")


def rsphere_triaxial(ell: Ellipsoid = None, method: str = "mean") -> float:
    """computes triaxial average of the semimajor and semiminor axes of the ellipsoid

    Parameters
    ----------
    ell : Ellipsoid, optional
          reference ellipsoid
    method: str, optional
        "mean" or "norm"

    Returns
    -------
    radius: float
        radius of sphere
    """

    if ell is None:
        ell = Ellipsoid()

    if method == "mean":
        return (2 * ell.semimajor_axis + ell.semiminor_axis) / 3
    elif method == "norm":
        return (ell.semimajor_axis ** 2 * ell.semiminor_axis) ** (1 / 3)
    else:
        raise Exception("pymap3d.rsphere.rsphere_triaxial: method must be mean or norm")


def rsphere_biaxial(ell: Ellipsoid = None, method: str = "mean") -> float:
    """computes biaxial average of the semimajor and semiminor axes of the ellipsoid

    Parameters
    ----------
    ell : Ellipsoid, optional
          reference ellipsoid
    method: str, optional
        "mean" or "norm"

    Returns
    -------
    radius: float
        radius of sphere
    """

    if ell is None:
        ell = Ellipsoid()

    if method == "mean":
        return (ell.semimajor_axis + ell.semiminor_axis) / 2
    elif method == "norm":
        return sqrt(ell.semimajor_axis * ell.semiminor_axis)
    else:
        raise Exception("pymap3d.rsphere.rsphere_biaxial: method must be mean or norm")
