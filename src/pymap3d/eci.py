""" transforms involving ECI earth-centered inertial """

from datetime import datetime
import typing
from numpy import array, sin, cos, column_stack, empty, atleast_1d

try:
    from astropy.coordinates import GCRS, ITRS, EarthLocation, CartesianRepresentation
    from astropy.time import Time
    import astropy.units as u
except ImportError:
    Time = None

from .sidereal import greenwichsrt, juliandate

try:
    from numpy.typing import ArrayLike
except ImportError:
    ArrayLike = typing.Any

__all__ = ["eci2ecef", "ecef2eci"]


def eci2ecef(
    x: ArrayLike, y: ArrayLike, z: ArrayLike, time: datetime, *, use_astropy: bool = True
) -> typing.Tuple[ArrayLike, ArrayLike, ArrayLike]:
    """
    Observer => Point  ECI  =>  ECEF

    J2000 frame

    Parameters
    ----------
    x : ArrayLike
        ECI x-location [meters]
    y : ArrayLike
        ECI y-location [meters]
    z : ArrayLike
        ECI z-location [meters]
    time : datetime.datetime
        time of obsevation (UTC)
    use_astropy: bool, optional
        use AstroPy (much more accurate)

    Results
    -------
    x_ecef : ArrayLike
        x ECEF coordinate
    y_ecef : ArrayLike
        y ECEF coordinate
    z_ecef : ArrayLike
        z ECEF coordinate
    """

    if use_astropy and Time is not None:
        gcrs = GCRS(CartesianRepresentation(x * u.m, y * u.m, z * u.m), obstime=time)
        itrs = gcrs.transform_to(ITRS(obstime=time))

        x_ecef = itrs.x.value
        y_ecef = itrs.y.value
        z_ecef = itrs.z.value
    else:
        x = atleast_1d(x)
        y = atleast_1d(y)
        z = atleast_1d(z)
        gst = atleast_1d(greenwichsrt(juliandate(time)))
        assert x.shape == y.shape == z.shape
        assert x.size == gst.size

        eci = column_stack((x.ravel(), y.ravel(), z.ravel()))
        ecef = empty((x.size, 3))
        for i in range(eci.shape[0]):
            ecef[i, :] = R3(gst[i]) @ eci[i, :].T

        x_ecef = ecef[:, 0].reshape(x.shape)
        y_ecef = ecef[:, 1].reshape(y.shape)
        z_ecef = ecef[:, 2].reshape(z.shape)

    return x_ecef, y_ecef, z_ecef


def ecef2eci(
    x: ArrayLike, y: ArrayLike, z: ArrayLike, time: datetime, *, use_astropy: bool = True
) -> typing.Tuple[ArrayLike, ArrayLike, ArrayLike]:
    """
    Point => Point   ECEF => ECI

    J2000 frame

    Parameters
    ----------

    x : ArrayLike
        target x ECEF coordinate
    y : ArrayLike
        target y ECEF coordinate
    z : ArrayLike
        target z ECEF coordinate
    time : datetime.datetime
        time of observation
    use_astropy: bool, optional
        use AstroPy (much more accurate)

    Results
    -------
    x_eci : ArrayLike
        x ECI coordinate
    y_eci : ArrayLike
        y ECI coordinate
    z_eci : ArrayLike
        z ECI coordinate
    """

    if use_astropy and Time is not None:
        itrs = ITRS(CartesianRepresentation(x * u.m, y * u.m, z * u.m), obstime=time)
        gcrs = itrs.transform_to(GCRS(obstime=time))
        eci = EarthLocation(*gcrs.cartesian.xyz)

        x_eci = eci.x.value
        y_eci = eci.y.value
        z_eci = eci.z.value
    else:
        x = atleast_1d(x)
        y = atleast_1d(y)
        z = atleast_1d(z)
        gst = atleast_1d(greenwichsrt(juliandate(time)))
        assert x.shape == y.shape == z.shape
        assert x.size == gst.size

        ecef = column_stack((x.ravel(), y.ravel(), z.ravel()))
        eci = empty((x.size, 3))
        for i in range(x.size):
            eci[i, :] = R3(gst[i]).T @ ecef[i, :]

        x_eci = eci[:, 0].reshape(x.shape)
        y_eci = eci[:, 1].reshape(y.shape)
        z_eci = eci[:, 2].reshape(z.shape)

    return x_eci, y_eci, z_eci


def R3(x: float) -> ArrayLike:
    """Rotation matrix for ECI"""
    return array([[cos(x), sin(x), 0], [-sin(x), cos(x), 0], [0, 0, 1]])
