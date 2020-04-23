""" transforms involving ECI earth-centered inertial """

from datetime import datetime
import typing

try:
    from astropy.coordinates import GCRS, ITRS, EarthLocation, CartesianRepresentation
    from astropy.time import Time
    import astropy.units as u
except ImportError:
    Time = None

__all__ = ["eci2ecef", "ecef2eci"]

if typing.TYPE_CHECKING:
    from numpy import ndarray


def eci2ecef(x: "ndarray", y: "ndarray", z: "ndarray", time: datetime) -> typing.Tuple["ndarray", "ndarray", "ndarray"]:
    """
    Observer => Point  ECI  =>  ECEF

    J2000 frame

    Parameters
    ----------
    x : "ndarray"
        ECI x-location [meters]
    y : "ndarray"
        ECI y-location [meters]
    z : "ndarray"
        ECI z-location [meters]
    time : datetime.datetime
        time of obsevation (UTC)

    Results
    -------
    x : "ndarray"
        target x ECEF coordinate
    y : "ndarray"
        target y ECEF coordinate
    z : "ndarray"
        target z ECEF coordinate
    """

    if Time is None:
        raise ImportError("pip install astropy")

    gcrs = GCRS(CartesianRepresentation(x * u.m, y * u.m, z * u.m), obstime=time)
    itrs = gcrs.transform_to(ITRS(obstime=time))

    return itrs.x.value, itrs.y.value, itrs.z.value


def ecef2eci(x: "ndarray", y: "ndarray", z: "ndarray", time: datetime) -> typing.Tuple["ndarray", "ndarray", "ndarray"]:
    """
    Point => Point   ECEF => ECI

    J2000 frame

    Parameters
    ----------

    x : "ndarray"
        target x ECEF coordinate
    y : "ndarray"
        target y ECEF coordinate
    z : "ndarray"
        target z ECEF coordinate
    time : datetime.datetime
        time of observation

    Results
    -------
    x : "ndarray"
        target x ECI coordinate
    y : "ndarray"
        target y ECI coordinate
    z : "ndarray"
        target z ECI coordinate
    """

    if Time is None:
        raise ImportError("pip install astropy")

    itrs = ITRS(CartesianRepresentation(x * u.m, y * u.m, z * u.m), obstime=time)
    gcrs = itrs.transform_to(GCRS(obstime=time))
    ecef = EarthLocation(*gcrs.cartesian.xyz)

    return ecef.x.value, ecef.y.value, ecef.z.value
