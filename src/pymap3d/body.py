"""Transforms involving an aerospace body frame."""

from __future__ import annotations

from .ecef import ecef2enu, enu2ecef
from .ellipsoid import Ellipsoid
from .enu import enu2geodetic, geodetic2enu
from .frames import (
    _body2ecef_rotation,
    _body2enu_rotation,
    _body2ned_rotation,
    _body2nwu_rotation,
    _body2sez_rotation,
    _ecef2body_rotation,
    _enu2body_rotation,
    _matvec3,
    _ned2body_rotation,
    _nwu2body_rotation,
    _sez2body_rotation,
)
from .ned import ecef2ned, ned2ecef, ned2geodetic, geodetic2ned
from .nwu import ecef2nwu, geodetic2nwu, nwu2ecef, nwu2geodetic
from .sez import ecef2sez, geodetic2sez, sez2ecef, sez2geodetic

__all__ = [
    "body_matrix",
    "ned2body_matrix",
    "body2ned_matrix",
    "enu2body_matrix",
    "body2enu_matrix",
    "sez2body_matrix",
    "body2sez_matrix",
    "nwu2body_matrix",
    "body2nwu_matrix",
    "ned2body",
    "body2ned",
    "enu2body",
    "body2enu",
    "sez2body",
    "body2sez",
    "nwu2body",
    "body2nwu",
    "ecef2body",
    "body2ecef",
    "ecef2bodyv",
    "body2ecefv",
    "geodetic2body",
    "body2geodetic",
]


def body_matrix(yaw, pitch, roll, sequence: str = "zyx", deg: bool = True):
    """
    Rotation matrix that maps NED coordinates into the body frame.

    The body frame uses aerospace axes: x forward, y right, z down.
    Positive yaw turns right, positive pitch is nose up, and positive
    roll is right wing down.
    """

    return _ned2body_rotation(yaw, pitch, roll, sequence=sequence, deg=deg)


def ned2body_matrix(yaw, pitch, roll, sequence: str = "zyx", deg: bool = True):
    """Rotation matrix that maps NED coordinates into the body frame."""

    return _ned2body_rotation(yaw, pitch, roll, sequence=sequence, deg=deg)


def body2ned_matrix(yaw, pitch, roll, sequence: str = "zyx", deg: bool = True):
    """Rotation matrix that maps body coordinates into NED."""

    return _body2ned_rotation(yaw, pitch, roll, sequence=sequence, deg=deg)


def enu2body_matrix(yaw, pitch, roll, sequence: str = "zyx", deg: bool = True):
    """Rotation matrix that maps ENU coordinates into the body frame."""

    return _enu2body_rotation(yaw, pitch, roll, sequence=sequence, deg=deg)


def body2enu_matrix(yaw, pitch, roll, sequence: str = "zyx", deg: bool = True):
    """Rotation matrix that maps body coordinates into ENU."""

    return _body2enu_rotation(yaw, pitch, roll, sequence=sequence, deg=deg)


def sez2body_matrix(yaw, pitch, roll, sequence: str = "zyx", deg: bool = True):
    """Rotation matrix that maps SEZ coordinates into the body frame."""

    return _sez2body_rotation(yaw, pitch, roll, sequence=sequence, deg=deg)


def body2sez_matrix(yaw, pitch, roll, sequence: str = "zyx", deg: bool = True):
    """Rotation matrix that maps body coordinates into SEZ."""

    return _body2sez_rotation(yaw, pitch, roll, sequence=sequence, deg=deg)


def nwu2body_matrix(yaw, pitch, roll, sequence: str = "zyx", deg: bool = True):
    """Rotation matrix that maps NWU coordinates into the body frame."""

    return _nwu2body_rotation(yaw, pitch, roll, sequence=sequence, deg=deg)


def body2nwu_matrix(yaw, pitch, roll, sequence: str = "zyx", deg: bool = True):
    """Rotation matrix that maps body coordinates into NWU."""

    return _body2nwu_rotation(yaw, pitch, roll, sequence=sequence, deg=deg)


def ned2body(n, e, d, yaw, pitch, roll, sequence: str = "zyx", deg: bool = True):
    """Convert NED coordinates into body coordinates."""

    return _matvec3(ned2body_matrix(yaw, pitch, roll, sequence=sequence, deg=deg), (n, e, d))


def body2ned(xb, yb, zb, yaw, pitch, roll, sequence: str = "zyx", deg: bool = True):
    """Convert body coordinates into NED."""

    return _matvec3(body2ned_matrix(yaw, pitch, roll, sequence=sequence, deg=deg), (xb, yb, zb))


def enu2body(e, n, u, yaw, pitch, roll, sequence: str = "zyx", deg: bool = True):
    """Convert ENU coordinates into body coordinates."""

    return _matvec3(enu2body_matrix(yaw, pitch, roll, sequence=sequence, deg=deg), (e, n, u))


def body2enu(xb, yb, zb, yaw, pitch, roll, sequence: str = "zyx", deg: bool = True):
    """Convert body coordinates into ENU."""

    return _matvec3(body2enu_matrix(yaw, pitch, roll, sequence=sequence, deg=deg), (xb, yb, zb))


def sez2body(s, e, z, yaw, pitch, roll, sequence: str = "zyx", deg: bool = True):
    """Convert SEZ coordinates into body coordinates."""

    return _matvec3(sez2body_matrix(yaw, pitch, roll, sequence=sequence, deg=deg), (s, e, z))


def body2sez(xb, yb, zb, yaw, pitch, roll, sequence: str = "zyx", deg: bool = True):
    """Convert body coordinates into SEZ."""

    return _matvec3(body2sez_matrix(yaw, pitch, roll, sequence=sequence, deg=deg), (xb, yb, zb))


def nwu2body(n, w, u, yaw, pitch, roll, sequence: str = "zyx", deg: bool = True):
    """Convert NWU coordinates into body coordinates."""

    return _matvec3(nwu2body_matrix(yaw, pitch, roll, sequence=sequence, deg=deg), (n, w, u))


def body2nwu(xb, yb, zb, yaw, pitch, roll, sequence: str = "zyx", deg: bool = True):
    """Convert body coordinates into NWU."""

    return _matvec3(body2nwu_matrix(yaw, pitch, roll, sequence=sequence, deg=deg), (xb, yb, zb))


def ecef2body(
    x,
    y,
    z,
    lat0,
    lon0,
    h0,
    yaw,
    pitch,
    roll,
    ell: Ellipsoid | None = None,
    frame: str = "ned",
    sequence: str = "zyx",
    deg: bool = True,
):
    """Convert a target ECEF point into body coordinates."""

    if frame == "ned":
        n, e, d = ecef2ned(x, y, z, lat0, lon0, h0, ell, deg=deg)
        return ned2body(n, e, d, yaw, pitch, roll, sequence=sequence, deg=deg)
    if frame == "enu":
        e, n, u = ecef2enu(x, y, z, lat0, lon0, h0, ell, deg=deg)
        return enu2body(e, n, u, yaw, pitch, roll, sequence=sequence, deg=deg)
    if frame == "sez":
        s, e, z = ecef2sez(x, y, z, lat0, lon0, h0, ell, deg=deg)
        return sez2body(s, e, z, yaw, pitch, roll, sequence=sequence, deg=deg)
    if frame == "nwu":
        n, w, u = ecef2nwu(x, y, z, lat0, lon0, h0, ell, deg=deg)
        return nwu2body(n, w, u, yaw, pitch, roll, sequence=sequence, deg=deg)

    raise ValueError("frame must be one of 'ned', 'enu', 'sez', or 'nwu'")


def body2ecef(
    xb,
    yb,
    zb,
    lat0,
    lon0,
    h0,
    yaw,
    pitch,
    roll,
    ell: Ellipsoid | None = None,
    frame: str = "ned",
    sequence: str = "zyx",
    deg: bool = True,
):
    """Convert body coordinates into a target ECEF point."""

    if frame == "ned":
        n, e, d = body2ned(xb, yb, zb, yaw, pitch, roll, sequence=sequence, deg=deg)
        return ned2ecef(n, e, d, lat0, lon0, h0, ell, deg=deg)
    if frame == "enu":
        e, n, u = body2enu(xb, yb, zb, yaw, pitch, roll, sequence=sequence, deg=deg)
        return enu2ecef(e, n, u, lat0, lon0, h0, ell, deg=deg)
    if frame == "sez":
        s, e, z = body2sez(xb, yb, zb, yaw, pitch, roll, sequence=sequence, deg=deg)
        return sez2ecef(s, e, z, lat0, lon0, h0, ell, deg=deg)
    if frame == "nwu":
        n, w, u = body2nwu(xb, yb, zb, yaw, pitch, roll, sequence=sequence, deg=deg)
        return nwu2ecef(n, w, u, lat0, lon0, h0, ell, deg=deg)

    raise ValueError("frame must be one of 'ned', 'enu', 'sez', or 'nwu'")


def ecef2bodyv(
    x,
    y,
    z,
    lat0,
    lon0,
    yaw,
    pitch,
    roll,
    frame: str = "ned",
    sequence: str = "zyx",
    deg: bool = True,
):
    """Convert an ECEF vector into body coordinates."""

    rotation = _ecef2body_rotation(
        lat0, lon0, yaw, pitch, roll, sequence=sequence, frame=frame, deg=deg
    )
    return _matvec3(rotation, (x, y, z))


def body2ecefv(
    xb,
    yb,
    zb,
    lat0,
    lon0,
    yaw,
    pitch,
    roll,
    frame: str = "ned",
    sequence: str = "zyx",
    deg: bool = True,
):
    """Convert a body-frame vector into ECEF."""

    rotation = _body2ecef_rotation(
        lat0, lon0, yaw, pitch, roll, sequence=sequence, frame=frame, deg=deg
    )
    return _matvec3(rotation, (xb, yb, zb))


def geodetic2body(
    lat,
    lon,
    h,
    lat0,
    lon0,
    h0,
    yaw,
    pitch,
    roll,
    ell: Ellipsoid | None = None,
    frame: str = "ned",
    sequence: str = "zyx",
    deg: bool = True,
):
    """Convert a target geodetic point into body coordinates."""

    if frame == "ned":
        n, e, d = geodetic2ned(lat, lon, h, lat0, lon0, h0, ell, deg=deg)
        return ned2body(n, e, d, yaw, pitch, roll, sequence=sequence, deg=deg)
    if frame == "enu":
        e, n, u = geodetic2enu(lat, lon, h, lat0, lon0, h0, ell, deg=deg)
        return enu2body(e, n, u, yaw, pitch, roll, sequence=sequence, deg=deg)
    if frame == "sez":
        s, e, z = geodetic2sez(lat, lon, h, lat0, lon0, h0, ell, deg=deg)
        return sez2body(s, e, z, yaw, pitch, roll, sequence=sequence, deg=deg)
    if frame == "nwu":
        n, w, u = geodetic2nwu(lat, lon, h, lat0, lon0, h0, ell, deg=deg)
        return nwu2body(n, w, u, yaw, pitch, roll, sequence=sequence, deg=deg)

    raise ValueError("frame must be one of 'ned', 'enu', 'sez', or 'nwu'")


def body2geodetic(
    xb,
    yb,
    zb,
    lat0,
    lon0,
    h0,
    yaw,
    pitch,
    roll,
    ell: Ellipsoid | None = None,
    frame: str = "ned",
    sequence: str = "zyx",
    deg: bool = True,
):
    """Convert body coordinates into a target geodetic point."""

    if frame == "ned":
        n, e, d = body2ned(xb, yb, zb, yaw, pitch, roll, sequence=sequence, deg=deg)
        return ned2geodetic(n, e, d, lat0, lon0, h0, ell, deg=deg)
    if frame == "enu":
        e, n, u = body2enu(xb, yb, zb, yaw, pitch, roll, sequence=sequence, deg=deg)
        return enu2geodetic(e, n, u, lat0, lon0, h0, ell, deg=deg)
    if frame == "sez":
        s, e, z = body2sez(xb, yb, zb, yaw, pitch, roll, sequence=sequence, deg=deg)
        return sez2geodetic(s, e, z, lat0, lon0, h0, ell, deg=deg)
    if frame == "nwu":
        n, w, u = body2nwu(xb, yb, zb, yaw, pitch, roll, sequence=sequence, deg=deg)
        return nwu2geodetic(n, w, u, lat0, lon0, h0, ell, deg=deg)

    raise ValueError("frame must be one of 'ned', 'enu', 'sez', or 'nwu'")
