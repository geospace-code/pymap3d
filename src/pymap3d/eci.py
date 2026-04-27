"""transforms involving ECI earth-centered inertial"""

from __future__ import annotations

from datetime import datetime
import sys
import logging

try:
    import numpy
    import astropy.units as u
    from astropy.coordinates import GCRS, ITRS, CartesianRepresentation, EarthLocation
except ImportError:
    pass

from .earth_orientation import eci_to_ecef_matrix, matvec3, transpose3

__all__ = ["eci2ecef", "ecef2eci"]


def eci2ecef(
    x,
    y,
    z,
    time: datetime,
    force_non_astropy: bool = False,
    *,
    delta_ut1: float = 0.0,
    xp: float = 0.0,
    yp: float = 0.0,
) -> tuple:
    """
    Observer => Point  ECI  =>  ECEF

    For the Astropy path, this converts from GCRS to ITRS.
    For the pure-Python path, this uses a J2000-like inertial frame with
    precession, truncated nutation, sidereal rotation, and optional polar motion.

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
    force_non_astropy : bool
        if True, force use of less accurate Numpy implementation even if Astropy is available
    delta_ut1 : float, optional
        UT1-UTC in seconds for the pure-Python path. Defaults to ``0.0``.
    xp : float, optional
        Polar motion x coordinate in arcseconds for the pure-Python path.
    yp : float, optional
        Polar motion y coordinate in arcseconds for the pure-Python path.

    Results
    -------
    x_ecef : float
        x ECEF coordinate
    y_ecef : float
        y ECEF coordinate
    z_ecef : float
        z ECEF coordinate
    """

    if "astropy" in sys.modules and not force_non_astropy:
        xe, ye, ze = eci2ecef_astropy(x, y, z, time)
    else:
        logging.warning(
            f"{__name__}: using pure-Python ECI/ECEF fallback, less accurate than Astropy"
        )
        xe, ye, ze = eci2ecef_numpy(x, y, z, time, delta_ut1=delta_ut1, xp=xp, yp=yp)

    return _squeeze_xyz(xe, ye, ze)


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


def eci2ecef_numpy(
    x,
    y,
    z,
    t: datetime,
    *,
    delta_ut1: float = 0.0,
    xp: float = 0.0,
    yp: float = 0.0,
) -> tuple:
    """
    eci2ecef using the internal pure-Python Earth orientation model.

    see eci2ecef() for description
    """

    rotation = eci_to_ecef_matrix(t, delta_ut1=delta_ut1, xp=xp, yp=yp)
    return _squeeze_xyz(*_apply_rotation(rotation, x, y, z))


def ecef2eci(
    x,
    y,
    z,
    time: datetime,
    force_non_astropy: bool = False,
    *,
    delta_ut1: float = 0.0,
    xp: float = 0.0,
    yp: float = 0.0,
) -> tuple:
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
    force_non_astropy : bool
        if True, force use of less accurate Numpy implementation even if Astropy is available
    delta_ut1 : float, optional
        UT1-UTC in seconds for the pure-Python path. Defaults to ``0.0``.
    xp : float, optional
        Polar motion x coordinate in arcseconds for the pure-Python path.
    yp : float, optional
        Polar motion y coordinate in arcseconds for the pure-Python path.

    Results
    -------
    x_eci : float
        x ECI coordinate
    y_eci : float
        y ECI coordinate
    z_eci : float
        z ECI coordinate
    """

    if "astropy" in sys.modules and not force_non_astropy:
        xe, ye, ze = ecef2eci_astropy(x, y, z, time)
    else:
        logging.warning(
            f"{__name__}: using pure-Python ECI/ECEF fallback, less accurate than Astropy"
        )
        xe, ye, ze = ecef2eci_numpy(x, y, z, time, delta_ut1=delta_ut1, xp=xp, yp=yp)

    return _squeeze_xyz(xe, ye, ze)


def ecef2eci_astropy(x, y, z, t: datetime) -> tuple:
    """ecef2eci using Astropy
    see ecef2eci() for description
    """
    itrs = ITRS(CartesianRepresentation(x * u.m, y * u.m, z * u.m), obstime=t)
    gcrs = itrs.transform_to(GCRS(obstime=t))
    eci = EarthLocation(*gcrs.cartesian.xyz)

    return eci.x.value, eci.y.value, eci.z.value


def ecef2eci_numpy(
    x,
    y,
    z,
    t: datetime,
    *,
    delta_ut1: float = 0.0,
    xp: float = 0.0,
    yp: float = 0.0,
) -> tuple:
    """ecef2eci using Numpy
    see ecef2eci() for description
    """

    rotation = transpose3(eci_to_ecef_matrix(t, delta_ut1=delta_ut1, xp=xp, yp=yp))
    return _squeeze_xyz(*_apply_rotation(rotation, x, y, z))


def _apply_rotation(rotation, x, y, z):
    """Apply a 3x3 rotation to scalar or NumPy coordinate inputs."""

    if "numpy" not in sys.modules:
        return matvec3(rotation, (x, y, z))

    x = numpy.atleast_1d(x)
    y = numpy.atleast_1d(y)
    z = numpy.atleast_1d(z)
    assert x.shape == y.shape == z.shape, (
        f"shape mismatch: x: {x.shape}  y: {y.shape}  z: {z.shape}"
    )

    vin = numpy.column_stack((x.ravel(), y.ravel(), z.ravel()))
    vout = numpy.empty((x.size, 3))
    for i in range(vin.shape[0]):
        vout[i, :] = matvec3(rotation, tuple(vin[i, :]))

    return (
        vout[:, 0].reshape(x.shape),
        vout[:, 1].reshape(y.shape),
        vout[:, 2].reshape(z.shape),
    )


def _squeeze_xyz(x, y, z):
    """Preserve scalar inputs while supporting array outputs."""

    try:
        return x.squeeze()[()], y.squeeze()[()], z.squeeze()[()]
    except AttributeError:
        return x, y, z
