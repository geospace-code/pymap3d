"""transforms involving ECI earth-centered inertial"""

from __future__ import annotations

from datetime import datetime
import sys
import logging

from ._typing import FloatLike, NDArray

try:
    import numpy as np
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

    if "astropy" in sys.modules and not force_non_astropy:
        xe, ye, ze = eci2ecef_astropy(x, y, z, time)
    elif "numpy" in sys.modules:
        logging.warning(f"{__name__}: Numpy implementation has much less accuracy than Astropy")
        xe, ye, ze = eci2ecef_numpy(x, y, z, time)
    else:
        raise ImportError("eci2ecef requires either Numpy or Astropy")

    return xe, ye, ze


def eci2ecef_astropy(
    x: FloatLike, y: FloatLike, z: FloatLike, t: datetime
) -> tuple[NDArray, NDArray, NDArray]:
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


def eci2ecef_numpy(x, y, z, t: datetime) -> tuple:
    """
    eci2ecef using Numpy

    see eci2ecef() for description
    """

    gst = greenwichsrt(juliandate(t))

    c = np.cos(gst)
    s = np.sin(gst)

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

    if "astropy" in sys.modules and not force_non_astropy:
        xe, ye, ze = ecef2eci_astropy(x, y, z, time)
    elif "numpy" in sys.modules:
        logging.warning(f"{__name__}: Numpy implementation has much less accuracy than Astropy")
        xe, ye, ze = ecef2eci_numpy(x, y, z, time)
    else:
        raise ImportError("ecef2eci requires either Numpy or Astropy")

    return xe, ye, ze


def ecef2eci_astropy(
    x: FloatLike, y: FloatLike, z: FloatLike, t: datetime
) -> tuple[NDArray, NDArray, NDArray]:
    """ecef2eci using Astropy
    see ecef2eci() for description
    """
    itrs = ITRS(CartesianRepresentation(x * u.m, y * u.m, z * u.m), obstime=t)
    gcrs = itrs.transform_to(GCRS(obstime=t))
    eci = EarthLocation(*gcrs.cartesian.xyz)

    return eci.x.value, eci.y.value, eci.z.value


def ecef2eci_numpy(x, y, z, t: datetime) -> tuple:
    """ecef2eci using Numpy
    see ecef2eci() for description
    """

    gst = greenwichsrt(juliandate(t))

    c = np.cos(gst)
    s = np.sin(gst)

    x_eci = c * x - s * y
    y_eci = s * x + c * y

    return x_eci, y_eci, z
