""" compute radii of auxiliary spheres"""
try:
    from numpy import radians, sin, cos, log, sqrt, degrees, asarray
except ImportError:
    from math import radians, sin, cos, log, sqrt, degrees
    asarray = None
from .ellipsoid import Ellipsoid
from .rcurve import rcurve_meridian, rcurve_transverse
from .vincenty import vdist

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
    """computes the radius of the sphere with equal volume as the ellipsoid"""
    if ell is None:
        ell = Ellipsoid()

    f = ell.flattening

    return ell.semimajor_axis * (1 - f / 3 - f ** 2 / 9)


def rsphere_authalic(ell: Ellipsoid = None) -> float:
    """computes the radius of the sphere with equal surface area as the ellipsoid"""
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
    """computes the radius of the sphere with equal meridional distances as the ellipsoid"""
    if ell is None:
        ell = Ellipsoid()
    return ((ell.semimajor_axis ** (3 / 2) + ell.semiminor_axis ** (3 / 2)) / 2) ** (2 / 3)


def rsphere_euler(lat1, lon1, lat2, lon2, ell: Ellipsoid = None, deg: bool = True) -> float:
    """computes the Euler radii of curvature at the midpoint of the
     great circle arc defined by the endpoints (lat1,lon1) and (lat2,lon2)"""
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


def rsphere_curve(lat, ell: Ellipsoid = None, deg: bool = True, method="mean") -> float:
    """computes the arithmetic average of the transverse and meridional
    radii of curvature at a specified latitude point"""

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


def rsphere_triaxial(ell: Ellipsoid = None, method="mean") -> float:
    """computes triaxial average of the semimajor and semiminor axes of the ellipsoid"""

    if ell is None:
        ell = Ellipsoid()

    if method == "mean":
        return (2 * ell.semimajor_axis + ell.semiminor_axis) / 3
    elif method == "norm":
        return (ell.semimajor_axis ** 2 * ell.semiminor_axis) ** (1 / 3)
    else:
        raise Exception("pymap3d.rsphere.rsphere_triaxial: method must be mean or norm")


def rsphere_biaxial(ell: Ellipsoid = None, method="mean") -> float:
    """computes biaxial average of the semimajor and semiminor axes of the ellipsoid"""

    if ell is None:
        ell = Ellipsoid()

    if method == "mean":
        return (ell.semimajor_axis + ell.semiminor_axis) / 2
    elif method == "norm":
        return sqrt(ell.semimajor_axis * ell.semiminor_axis)
    else:
        raise Exception("pymap3d.rsphere.rsphere_biaxial: method must be mean or norm")
