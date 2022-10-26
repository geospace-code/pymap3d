""" compute radii of auxiliary spheres"""

from __future__ import annotations

try:
    from numpy import asarray
except ImportError:
    pass

from . import rcurve
from .ellipsoid import Ellipsoid
from .mathfun import cos, degrees, log, radians, sin, sqrt
from .vincenty import vdist

__all__ = [
    "eqavol",
    "authalic",
    "rectifying",
    "euler",
    "curve",
    "triaxial",
    "biaxial",
]


def eqavol(ell: Ellipsoid = None) -> float:
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
        ell = Ellipsoid.from_name("wgs84")

    f = ell.flattening

    return ell.semimajor_axis * (1 - f / 3 - f**2 / 9)


def authalic(ell: Ellipsoid = None) -> float:
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
        ell = Ellipsoid.from_name("wgs84")

    e = ell.eccentricity

    if e > 0:
        f1 = ell.semimajor_axis**2 / 2
        f2 = (1 - e**2) / (2 * e)
        f3 = log((1 + e) / (1 - e))
        return sqrt(f1 * (1 + f2 * f3))
    else:
        return ell.semimajor_axis


def rectifying(ell: Ellipsoid = None) -> float:
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
        ell = Ellipsoid.from_name("wgs84")
    return ((ell.semimajor_axis ** (3 / 2) + ell.semiminor_axis ** (3 / 2)) / 2) ** (2 / 3)


def euler(
    lat1,
    lon1,
    lat2,
    lon2,
    ell: Ellipsoid = None,
    deg: bool = True,
):
    """computes the Euler radii of curvature at the midpoint of the
     great circle arc defined by the endpoints (lat1,lon1) and (lat2,lon2)

    Parameters
    ----------
    lat1, lat2 : float
        geodetic latitudes (degrees)
    lon1, lon2 : float
       geodetic longitudes (degrees)
    ell : Ellipsoid, optional
          reference ellipsoid
    deg : bool, optional
          degrees input/output  (False: radians in/out)

    Returns
    -------
    radius: float
        radius of sphere
    """
    if not deg:
        lat1, lon1, lat2, lon2 = degrees(lat1), degrees(lon1), degrees(lat2), degrees(lon2)

    try:
        lat1, lat2 = asarray(lat1), asarray(lat2)
    except NameError:
        pass

    latmid = lat1 + (lat2 - lat1) / 2  # compute the midpoint

    # compute azimuth
    az = vdist(lat1, lon1, lat2, lon2, ell=ell)[1]

    #   compute meridional and transverse radii of curvature
    rho = rcurve.meridian(latmid, ell, deg=True)
    nu = rcurve.transverse(latmid, ell, deg=True)

    az = radians(az)
    den = rho * sin(az) ** 2 + nu * cos(az) ** 2

    #  compute radius of the arc from point 1 to point 2
    return rho * nu / den


def curve(lat, ell: Ellipsoid = None, deg: bool = True, method: str = "mean"):
    """computes the arithmetic average of the transverse and meridional
    radii of curvature at a specified latitude point

    Parameters
    ----------
    lat1 : float
        geodetic latitudes (degrees)
    ell : Ellipsoid, optional
          reference ellipsoid
    method: str, optional
        "mean" or "norm"
    deg : bool, optional
          degrees input/output  (False: radians in/out)

    Returns
    -------
    radius: float
        radius of sphere
    """

    if deg:
        lat = radians(lat)

    rho = rcurve.meridian(lat, ell, deg=False)
    nu = rcurve.transverse(lat, ell, deg=False)

    if method == "mean":
        return (rho + nu) / 2
    elif method == "norm":
        return sqrt(rho * nu)
    else:
        raise ValueError("method must be mean or norm")


def triaxial(ell: Ellipsoid = None, method: str = "mean") -> float:
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
        ell = Ellipsoid.from_name("wgs84")

    if method == "mean":
        return (2 * ell.semimajor_axis + ell.semiminor_axis) / 3
    elif method == "norm":
        return (ell.semimajor_axis**2 * ell.semiminor_axis) ** (1 / 3)
    else:
        raise ValueError("method must be mean or norm")


def biaxial(ell: Ellipsoid = None, method: str = "mean") -> float:
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
        ell = Ellipsoid.from_name("wgs84")

    if method == "mean":
        return (ell.semimajor_axis + ell.semiminor_axis) / 2
    elif method == "norm":
        return sqrt(ell.semimajor_axis * ell.semiminor_axis)
    else:
        raise ValueError("method must be mean or norm")
