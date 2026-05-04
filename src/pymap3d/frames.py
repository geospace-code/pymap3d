"""Shared local-frame rotation helpers."""

from __future__ import annotations

from .mathfun import cos, radians, sin


def _matmul3(a, b):
    """Multiply two 3x3 matrices represented as tuples."""

    return tuple(
        tuple(sum(a[i][k] * b[k][j] for k in range(3)) for j in range(3)) for i in range(3)
    )


def _transpose3(matrix):
    """Transpose a 3x3 matrix represented as tuples."""

    return tuple(tuple(matrix[j][i] for j in range(3)) for i in range(3))


def _matvec3(matrix, vector):
    """Apply a 3x3 matrix to a 3-vector."""

    return tuple(sum(matrix[i][j] * vector[j] for j in range(3)) for i in range(3))


def _rotation_x(angle):
    """Active rotation matrix about x."""

    cos_angle = cos(angle)
    sin_angle = sin(angle)
    return (
        (1.0, 0.0, 0.0),
        (0.0, cos_angle, -sin_angle),
        (0.0, sin_angle, cos_angle),
    )


def _rotation_y(angle):
    """Active rotation matrix about y."""

    cos_angle = cos(angle)
    sin_angle = sin(angle)
    return (
        (cos_angle, 0.0, sin_angle),
        (0.0, 1.0, 0.0),
        (-sin_angle, 0.0, cos_angle),
    )


def _rotation_z(angle):
    """Active rotation matrix about z."""

    cos_angle = cos(angle)
    sin_angle = sin(angle)
    return (
        (cos_angle, -sin_angle, 0.0),
        (sin_angle, cos_angle, 0.0),
        (0.0, 0.0, 1.0),
    )


def _ecef2enu_rotation(lat0, lon0, deg: bool = True):
    """Rotation matrix that maps ECEF vectors into ENU coordinates."""

    if deg:
        lat0 = radians(lat0)
        lon0 = radians(lon0)

    sin_lat = sin(lat0)
    cos_lat = cos(lat0)
    sin_lon = sin(lon0)
    cos_lon = cos(lon0)

    return (
        (-sin_lon, cos_lon, 0.0),
        (-sin_lat * cos_lon, -sin_lat * sin_lon, cos_lat),
        (cos_lat * cos_lon, cos_lat * sin_lon, sin_lat),
    )


def _enu2ecef_rotation(lat0, lon0, deg: bool = True):
    """Rotation matrix that maps ENU vectors into ECEF coordinates."""

    return _transpose3(_ecef2enu_rotation(lat0, lon0, deg=deg))


def _ecef2ned_rotation(lat0, lon0, deg: bool = True):
    """Rotation matrix that maps ECEF vectors into NED coordinates."""

    return _matmul3(_ENU2NED, _ecef2enu_rotation(lat0, lon0, deg=deg))


def _ned2ecef_rotation(lat0, lon0, deg: bool = True):
    """Rotation matrix that maps NED vectors into ECEF coordinates."""

    return _transpose3(_ecef2ned_rotation(lat0, lon0, deg=deg))


def _ecef2sez_rotation(lat0, lon0, deg: bool = True):
    """Rotation matrix that maps ECEF vectors into SEZ coordinates."""

    return _matmul3(_ENU2SEZ, _ecef2enu_rotation(lat0, lon0, deg=deg))


def _sez2ecef_rotation(lat0, lon0, deg: bool = True):
    """Rotation matrix that maps SEZ vectors into ECEF coordinates."""

    return _transpose3(_ecef2sez_rotation(lat0, lon0, deg=deg))


def _ecef2nwu_rotation(lat0, lon0, deg: bool = True):
    """Rotation matrix that maps ECEF vectors into NWU coordinates."""

    return _matmul3(_ENU2NWU, _ecef2enu_rotation(lat0, lon0, deg=deg))


def _nwu2ecef_rotation(lat0, lon0, deg: bool = True):
    """Rotation matrix that maps NWU vectors into ECEF coordinates."""

    return _transpose3(_ecef2nwu_rotation(lat0, lon0, deg=deg))


def _ned2body_rotation(yaw, pitch, roll, sequence: str = "zyx", deg: bool = True):
    """
    Rotation matrix that maps NED vectors into body coordinates.

    The body frame uses aerospace axes: x forward, y right, z down.
    Positive yaw turns right, positive pitch is nose up, and positive
    roll is right wing down.
    """

    if sequence != "zyx":
        raise ValueError("only sequence='zyx' is currently supported")

    if deg:
        yaw = radians(yaw)
        pitch = radians(pitch)
        roll = radians(roll)

    body2ned = _matmul3(_rotation_z(yaw), _matmul3(_rotation_y(pitch), _rotation_x(roll)))
    return _transpose3(body2ned)


def _body2ned_rotation(yaw, pitch, roll, sequence: str = "zyx", deg: bool = True):
    """Rotation matrix that maps body vectors into NED coordinates."""

    return _transpose3(_ned2body_rotation(yaw, pitch, roll, sequence=sequence, deg=deg))


def _enu2body_rotation(yaw, pitch, roll, sequence: str = "zyx", deg: bool = True):
    """Rotation matrix that maps ENU vectors into body coordinates."""

    return _matmul3(_ned2body_rotation(yaw, pitch, roll, sequence=sequence, deg=deg), _ENU2NED)


def _body2enu_rotation(yaw, pitch, roll, sequence: str = "zyx", deg: bool = True):
    """Rotation matrix that maps body vectors into ENU coordinates."""

    return _transpose3(_enu2body_rotation(yaw, pitch, roll, sequence=sequence, deg=deg))


def _sez2body_rotation(yaw, pitch, roll, sequence: str = "zyx", deg: bool = True):
    """Rotation matrix that maps SEZ vectors into body coordinates."""

    return _matmul3(_enu2body_rotation(yaw, pitch, roll, sequence=sequence, deg=deg), _transpose3(_ENU2SEZ))


def _body2sez_rotation(yaw, pitch, roll, sequence: str = "zyx", deg: bool = True):
    """Rotation matrix that maps body vectors into SEZ coordinates."""

    return _transpose3(_sez2body_rotation(yaw, pitch, roll, sequence=sequence, deg=deg))


def _nwu2body_rotation(yaw, pitch, roll, sequence: str = "zyx", deg: bool = True):
    """Rotation matrix that maps NWU vectors into body coordinates."""

    return _matmul3(_enu2body_rotation(yaw, pitch, roll, sequence=sequence, deg=deg), _transpose3(_ENU2NWU))


def _body2nwu_rotation(yaw, pitch, roll, sequence: str = "zyx", deg: bool = True):
    """Rotation matrix that maps body vectors into NWU coordinates."""

    return _transpose3(_nwu2body_rotation(yaw, pitch, roll, sequence=sequence, deg=deg))


def _ecef2body_rotation(
    lat0,
    lon0,
    yaw,
    pitch,
    roll,
    sequence: str = "zyx",
    frame: str = "ned",
    deg: bool = True,
):
    """Rotation matrix that maps ECEF vectors into body coordinates."""

    if frame == "ned":
        local2body = _ned2body_rotation(yaw, pitch, roll, sequence=sequence, deg=deg)
        ecef2local = _ecef2ned_rotation(lat0, lon0, deg=deg)
    elif frame == "enu":
        local2body = _enu2body_rotation(yaw, pitch, roll, sequence=sequence, deg=deg)
        ecef2local = _ecef2enu_rotation(lat0, lon0, deg=deg)
    elif frame == "sez":
        local2body = _sez2body_rotation(yaw, pitch, roll, sequence=sequence, deg=deg)
        ecef2local = _ecef2sez_rotation(lat0, lon0, deg=deg)
    elif frame == "nwu":
        local2body = _nwu2body_rotation(yaw, pitch, roll, sequence=sequence, deg=deg)
        ecef2local = _ecef2nwu_rotation(lat0, lon0, deg=deg)
    else:
        raise ValueError("frame must be one of 'ned', 'enu', 'sez', or 'nwu'")

    return _matmul3(local2body, ecef2local)


def _body2ecef_rotation(
    lat0,
    lon0,
    yaw,
    pitch,
    roll,
    sequence: str = "zyx",
    frame: str = "ned",
    deg: bool = True,
):
    """Rotation matrix that maps body vectors into ECEF coordinates."""

    return _transpose3(
        _ecef2body_rotation(
            lat0,
            lon0,
            yaw,
            pitch,
            roll,
            sequence=sequence,
            frame=frame,
            deg=deg,
        )
    )


_ENU2NED = (
    (0.0, 1.0, 0.0),
    (1.0, 0.0, 0.0),
    (0.0, 0.0, -1.0),
)

_ENU2SEZ = (
    (0.0, -1.0, 0.0),
    (1.0, 0.0, 0.0),
    (0.0, 0.0, 1.0),
)

_ENU2NWU = (
    (0.0, 1.0, 0.0),
    (-1.0, 0.0, 0.0),
    (0.0, 0.0, 1.0),
)
