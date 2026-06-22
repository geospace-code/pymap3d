"""transforms involving ECI earth-centered inertial"""

from __future__ import annotations

from datetime import datetime

import sys

try:
    import numpy
    import astropy.units as u
    from astropy.coordinates import (
        GCRS,
        ITRS,
        CartesianDifferential,
        CartesianRepresentation,
        EarthLocation,
    )
except ImportError:
    pass

from .earth_orientation import eci_to_ecef_matrix, eci_to_ecef_matrix_rate, matvec3, transpose3

__all__ = ["eci2ecef", "ecef2eci", "eci2ecef_state", "ecef2eci_state"]


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
    delta_ut1 : float, optional
        UT1-UTC in seconds for the pure-Python path. Defaults to ``0.0``.
    xp : float, optional
        Polar motion x coordinate in arcseconds for the pure-Python path.
    yp : float, optional
        Polar motion y coordinate in arcseconds for the pure-Python path.

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
    else:
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
    delta_ut1 : float, optional
        UT1-UTC in seconds for the pure-Python path. Defaults to ``0.0``.
    xp : float, optional
        Polar motion x coordinate in arcseconds for the pure-Python path.
    yp : float, optional
        Polar motion y coordinate in arcseconds for the pure-Python path.

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
    else:
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


def eci2ecef_state(
    x,
    y,
    z,
    vx,
    vy,
    vz,
    time: datetime,
    force_non_astropy: bool = False,
    *,
    delta_ut1: float = 0.0,
    xp: float = 0.0,
    yp: float = 0.0,
) -> tuple:
    """
    Transform an ECI position/velocity state to ECEF.

    Velocity transforms depend on both position and velocity because ECEF rotates
    with the Earth.
    """

    if "astropy" in sys.modules and not force_non_astropy:
        x, y, z, vx, vy, vz = eci2ecef_state_astropy(x, y, z, vx, vy, vz, time)
    else:
        x, y, z, vx, vy, vz = eci2ecef_state_numpy(
            x,
            y,
            z,
            vx,
            vy,
            vz,
            time,
            delta_ut1=delta_ut1,
            xp=xp,
            yp=yp,
        )

    return _squeeze_state(x, y, z, vx, vy, vz)


def eci2ecef_state_astropy(x, y, z, vx, vy, vz, t: datetime) -> tuple:
    """eci2ecef_state using Astropy."""

    gcrs = GCRS(
        CartesianRepresentation(
            x * u.m,
            y * u.m,
            z * u.m,
            differentials=CartesianDifferential(vx * u.m / u.s, vy * u.m / u.s, vz * u.m / u.s),
        ),
        obstime=t,
    )
    itrs = gcrs.transform_to(ITRS(obstime=t))
    rep = itrs.cartesian
    diff = rep.differentials["s"]

    return (
        rep.x.to_value(u.m),
        rep.y.to_value(u.m),
        rep.z.to_value(u.m),
        diff.d_x.to_value(u.m / u.s),
        diff.d_y.to_value(u.m / u.s),
        diff.d_z.to_value(u.m / u.s),
    )


def eci2ecef_state_numpy(
    x,
    y,
    z,
    vx,
    vy,
    vz,
    t: datetime,
    *,
    delta_ut1: float = 0.0,
    xp: float = 0.0,
    yp: float = 0.0,
) -> tuple:
    """eci2ecef_state using the internal pure-Python Earth orientation model."""

    rotation = eci_to_ecef_matrix(t, delta_ut1=delta_ut1, xp=xp, yp=yp)
    rotation_rate = eci_to_ecef_matrix_rate(t, delta_ut1=delta_ut1, xp=xp, yp=yp)
    return _apply_state_rotation(
            rotation,
            rotation_rate,
            x,
            y,
            z,
            vx,
            vy,
            vz,
            inverse=False,
        )


def ecef2eci_state(
    x,
    y,
    z,
    vx,
    vy,
    vz,
    time: datetime,
    force_non_astropy: bool = False,
    *,
    delta_ut1: float = 0.0,
    xp: float = 0.0,
    yp: float = 0.0,
) -> tuple:
    """
    Transform an ECEF position/velocity state to ECI.

    Velocity transforms depend on both position and velocity because ECEF rotates
    with the Earth.
    """

    if "astropy" in sys.modules and not force_non_astropy:
        x, y, z, vx, vy, vz = ecef2eci_state_astropy(x, y, z, vx, vy, vz, time)
    else:
        x, y, z, vx, vy, vz = ecef2eci_state_numpy(
            x,
            y,
            z,
            vx,
            vy,
            vz,
            time,
            delta_ut1=delta_ut1,
            xp=xp,
            yp=yp,
        )

    return _squeeze_state(x, y, z, vx, vy, vz)


def ecef2eci_state_astropy(x, y, z, vx, vy, vz, t: datetime) -> tuple:
    """ecef2eci_state using Astropy."""

    itrs = ITRS(
        CartesianRepresentation(
            x * u.m,
            y * u.m,
            z * u.m,
            differentials=CartesianDifferential(vx * u.m / u.s, vy * u.m / u.s, vz * u.m / u.s),
        ),
        obstime=t,
    )
    gcrs = itrs.transform_to(GCRS(obstime=t))
    rep = gcrs.cartesian
    diff = rep.differentials["s"]

    return (
        rep.x.to_value(u.m),
        rep.y.to_value(u.m),
        rep.z.to_value(u.m),
        diff.d_x.to_value(u.m / u.s),
        diff.d_y.to_value(u.m / u.s),
        diff.d_z.to_value(u.m / u.s),
    )


def ecef2eci_state_numpy(
    x,
    y,
    z,
    vx,
    vy,
    vz,
    t: datetime,
    *,
    delta_ut1: float = 0.0,
    xp: float = 0.0,
    yp: float = 0.0,
) -> tuple:
    """ecef2eci_state using the internal pure-Python Earth orientation model."""

    rotation = eci_to_ecef_matrix(t, delta_ut1=delta_ut1, xp=xp, yp=yp)
    rotation_rate = eci_to_ecef_matrix_rate(t, delta_ut1=delta_ut1, xp=xp, yp=yp)
    return _apply_state_rotation(
            rotation,
            rotation_rate,
            x,
            y,
            z,
            vx,
            vy,
            vz,
            inverse=True,
        )


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


def _apply_state_rotation(
    rotation,
    rotation_rate,
    x,
    y,
    z,
    vx,
    vy,
    vz,
    *,
    inverse: bool,
):
    """Apply a position/velocity frame transform."""

    if "numpy" not in sys.modules:
        return _apply_state_rotation_scalar(
            rotation, rotation_rate, x, y, z, vx, vy, vz, inverse=inverse
        )

    x = numpy.atleast_1d(x)
    y = numpy.atleast_1d(y)
    z = numpy.atleast_1d(z)
    vx = numpy.atleast_1d(vx)
    vy = numpy.atleast_1d(vy)
    vz = numpy.atleast_1d(vz)
    assert x.shape == y.shape == z.shape == vx.shape == vy.shape == vz.shape, (
        "shape mismatch between position and velocity inputs"
    )

    pos_in = numpy.column_stack((x.ravel(), y.ravel(), z.ravel()))
    vel_in = numpy.column_stack((vx.ravel(), vy.ravel(), vz.ravel()))
    pos_out = numpy.empty((x.size, 3))
    vel_out = numpy.empty((x.size, 3))

    for i in range(pos_in.shape[0]):
        pos_i, vel_i = _apply_state_rotation_scalar(
            rotation,
            rotation_rate,
            *tuple(pos_in[i, :]),
            *tuple(vel_in[i, :]),
            inverse=inverse,
        )
        pos_out[i, :] = pos_i
        vel_out[i, :] = vel_i

    return (
        pos_out[:, 0].reshape(x.shape),
        pos_out[:, 1].reshape(y.shape),
        pos_out[:, 2].reshape(z.shape),
        vel_out[:, 0].reshape(vx.shape),
        vel_out[:, 1].reshape(vy.shape),
        vel_out[:, 2].reshape(vz.shape),
    )


def _apply_state_rotation_scalar(
    rotation,
    rotation_rate,
    x,
    y,
    z,
    vx,
    vy,
    vz,
    *,
    inverse: bool,
):
    """Scalar state rotation helper."""

    if inverse:
        position = matvec3(transpose3(rotation), (x, y, z))
        coriolis = matvec3(rotation_rate, position)
        velocity = matvec3(
            transpose3(rotation),
            (vx - coriolis[0], vy - coriolis[1], vz - coriolis[2]),
        )
    else:
        position = matvec3(rotation, (x, y, z))
        spin = matvec3(rotation_rate, (x, y, z))
        rotated_velocity = matvec3(rotation, (vx, vy, vz))
        velocity = tuple(rotated_velocity[i] + spin[i] for i in range(3))

    return *position, *velocity


def _squeeze_xyz(x, y, z):
    """Preserve scalar inputs while supporting array outputs."""

    try:
        return x.squeeze()[()], y.squeeze()[()], z.squeeze()[()]
    except AttributeError:
        return x, y, z


def _squeeze_state(x, y, z, vx, vy, vz):
    """Preserve scalar state inputs while supporting array outputs."""

    try:
        return (
            x.squeeze()[()],
            y.squeeze()[()],
            z.squeeze()[()],
            vx.squeeze()[()],
            vy.squeeze()[()],
            vz.squeeze()[()],
        )
    except AttributeError:
        return x, y, z, vx, vy, vz
