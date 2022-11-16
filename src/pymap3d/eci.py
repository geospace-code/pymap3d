""" transforms involving ECI earth-centered inertial """

from __future__ import annotations

from datetime import datetime
from typing import Any, overload

from numpy import array, atleast_1d, column_stack, cos, empty, sin
from numpy.typing import NDArray

try:
    import astropy.units as u
    from astropy.coordinates import GCRS, ITRS, CartesianRepresentation, EarthLocation
except ImportError:
    pass
from ._types import ArrayLike
from .sidereal import greenwichsrt, juliandate

__all__ = ["eci2ecef", "ecef2eci"]


@overload
def eci2ecef(x: float, y: float, z: float, time: datetime) -> tuple[float, float, float]:
    pass


@overload
def eci2ecef(
    x: ArrayLike, y: ArrayLike, z: ArrayLike, time: datetime
) -> tuple[NDArray[Any], NDArray[Any], NDArray[Any]]:
    pass


def eci2ecef(
    x: float | ArrayLike,
    y: float | ArrayLike,
    z: float | ArrayLike,
    time: datetime,
) -> tuple[float, float, float] | tuple[NDArray[Any], NDArray[Any], NDArray[Any]]:
    """
    Observer => Point  ECI  =>  ECEF

    J2000 frame

    Parameters
    ----------
    x : float
        ECI x-location [meters]
    y : float
        ECI y-location [meters]
    z : float
        ECI z-location [meters]
    time : datetime.datetime
        time of obsevation (UTC)

    Results
    -------
    x_ecef : float
        x ECEF coordinate
    y_ecef : float
        y ECEF coordinate
    z_ecef : float
        z ECEF coordinate
    """

    try:
        gcrs = GCRS(CartesianRepresentation(x * u.m, y * u.m, z * u.m), obstime=time)
        itrs = gcrs.transform_to(ITRS(obstime=time))

        x_ecef = itrs.x.value
        y_ecef = itrs.y.value
        z_ecef = itrs.z.value
    except NameError:
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


@overload
def ecef2eci(x: float, y: float, z: float, time: datetime) -> tuple[float, float, float]:
    pass


@overload
def ecef2eci(
    x: ArrayLike, y: ArrayLike, z: ArrayLike, time: datetime
) -> tuple[NDArray[Any], NDArray[Any], NDArray[Any]]:
    pass


def ecef2eci(
    x: float | ArrayLike,
    y: float | ArrayLike,
    z: float | ArrayLike,
    time: datetime,
) -> tuple[float, float, float] | tuple[NDArray[Any], NDArray[Any], NDArray[Any]]:
    """
    Point => Point   ECEF => ECI

    J2000 frame

    Parameters
    ----------

    x : float
        target x ECEF coordinate
    y : float
        target y ECEF coordinate
    z : float
        target z ECEF coordinate
    time : datetime.datetime
        time of observation

    Results
    -------
    x_eci : float
        x ECI coordinate
    y_eci : float
        y ECI coordinate
    z_eci : float
        z ECI coordinate
    """

    try:
        itrs = ITRS(CartesianRepresentation(x * u.m, y * u.m, z * u.m), obstime=time)
        gcrs = itrs.transform_to(GCRS(obstime=time))
        eci = EarthLocation(*gcrs.cartesian.xyz)

        x_eci = eci.x.value
        y_eci = eci.y.value
        z_eci = eci.z.value
    except NameError:
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


def R3(x: float) -> NDArray[Any]:
    """Rotation matrix for ECI"""
    return array([[cos(x), sin(x), 0], [-sin(x), cos(x), 0], [0, 0, 1]])
