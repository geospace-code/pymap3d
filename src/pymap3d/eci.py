"""transforms involving ECI earth-centered inertial"""

from __future__ import annotations

from datetime import datetime

from .mathfun import cos, sin
from ._typing import FloatLike

try:
    import astropy.units as u
    from astropy.coordinates import GCRS, ITRS, CartesianRepresentation, EarthLocation
except ImportError:
    pass


from .sidereal import greenwichsrt, juliandate

__all__ = ["eci2ecef", "ecef2eci"]


def eci2ecef(
    x: FloatLike, y: FloatLike, z: FloatLike, time: datetime, force_non_astropy: bool = False
) -> tuple:
    """
    Observer => Point  ECI  =>  ECEF

    J2000 frame

    Parameters
    ----------
    x : array-like float
        ECI x-location [meters]
    y : array-like float
        ECI y-location [meters]
    z : array-like float
        ECI z-location [meters]
    time : datetime.datetime
        time of obsevation (UTC)
    force_non_astropy : bool
        if True, force use of less accurate Numpy implementation even if Astropy is available

    Results
    -------
    x_ecef : array-like float
        x ECEF coordinate
    y_ecef : array-like float
        y ECEF coordinate
    z_ecef : array-like float
        z ECEF coordinate
    """

    if force_non_astropy:
        return eci2ecef_stdlib(x, y, z, time)

    try:
        return eci2ecef_astropy(x, y, z, time)
    except NameError:
        return eci2ecef_stdlib(x, y, z, time)


def eci2ecef_astropy(x, y, z, t: datetime) -> tuple:
    """
    eci2ecef using Astropy

    see eci2ecef() for description
    """

    gcrs = GCRS(CartesianRepresentation(x * u.m, y * u.m, z * u.m), obstime=t)
    itrs = gcrs.transform_to(ITRS(obstime=t))

    x_ecef = itrs.x.value
    y_ecef = itrs.y.value
    z_ecef = itrs.z.value

    return x_ecef, y_ecef, z_ecef


def eci2ecef_stdlib(x, y, z, t: datetime) -> tuple:
    """ eci2ecef without Astropy
    see eci2ecef() for description
    """

    gst = greenwichsrt(juliandate(t))

    c = cos(gst)
    s = sin(gst)

    x_ecef = c * x + s * y
    y_ecef = -s * x + c * y

    return x_ecef, y_ecef, z


def ecef2eci(
    x: FloatLike, y: FloatLike, z: FloatLike, time: datetime, force_non_astropy: bool = False
) -> tuple:
    """
    Point => Point   ECEF => ECI

    J2000 frame

    Parameters
    ----------

    x : array-like float
        point x ECEF coordinate
    y : array-like float
        point y ECEF coordinate
    z : array-like float
        point z ECEF coordinate
    time : datetime.datetime
        time of observation
    force_non_astropy : bool
        if True, force use of less accurate Numpy implementation even if Astropy is available

    Results
    -------
    x_eci : array-like float
        x ECI coordinate
    y_eci : array-like float
        y ECI coordinate
    z_eci : array-like float
        z ECI coordinate
    """

    if force_non_astropy:
        return ecef2eci_stdlib(x, y, z, time)

    try:
        return ecef2eci_astropy(x, y, z, time)
    except NameError:
        return ecef2eci_stdlib(x, y, z, time)


def ecef2eci_astropy(x, y, z, t: datetime) -> tuple:
    """ecef2eci using Astropy
    see ecef2eci() for description
    """
    itrs = ITRS(CartesianRepresentation(x * u.m, y * u.m, z * u.m), obstime=t)
    gcrs = itrs.transform_to(GCRS(obstime=t))
    eci = EarthLocation(*gcrs.cartesian.xyz)

    return eci.x.value, eci.y.value, eci.z.value


def ecef2eci_stdlib(x, y, z, t: datetime) -> tuple:
    """ecef2eci without Astropy
    see ecef2eci() for description
    """

    gst = greenwichsrt(juliandate(t))

    c = cos(gst)
    s = sin(gst)

    x_eci = c * x - s * y
    y_eci = s * x + c * y

    return x_eci, y_eci, z
