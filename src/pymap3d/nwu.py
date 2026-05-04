"""Transforms involving NWU North West Up."""

from __future__ import annotations

from .ecef import ecef2enu, ecef2enuv, enu2ecef
from .ellipsoid import Ellipsoid
from .enu import aer2enu, enu2aer, enu2ecefv, geodetic2enu, enu2geodetic
from .frames import _ecef2nwu_rotation, _nwu2ecef_rotation

__all__ = [
    "aer2nwu",
    "nwu2aer",
    "nwu2geodetic",
    "nwu2ecef",
    "nwu2ecefv",
    "ecef2nwu",
    "ecef2nwuv",
    "geodetic2nwu",
    "ecef2nwu_matrix",
    "nwu2ecef_matrix",
]


def aer2nwu(az, elev, slantRange, deg: bool = True) -> tuple:
    """
    Convert azimuth, elevation, range to North, West, Up.
    """

    e, n, u = aer2enu(az, elev, slantRange, deg=deg)
    return n, -e, u


def nwu2aer(n, w, u, deg: bool = True) -> tuple:
    """
    Convert North, West, Up to azimuth, elevation, range.
    """

    return enu2aer(-w, n, u, deg=deg)


def nwu2geodetic(
    n,
    w,
    u,
    lat0,
    lon0,
    h0,
    ell: Ellipsoid | None = None,
    deg: bool = True,
) -> tuple:
    """
    Convert North, West, Up to target geodetic coordinates.
    """

    return enu2geodetic(-w, n, u, lat0, lon0, h0, ell, deg=deg)


def nwu2ecef(
    n,
    w,
    u,
    lat0,
    lon0,
    h0,
    ell: Ellipsoid | None = None,
    deg: bool = True,
) -> tuple:
    """
    Convert North, West, Up to target ECEF coordinates.
    """

    return enu2ecef(-w, n, u, lat0, lon0, h0, ell, deg=deg)


def nwu2ecefv(n, w, u, lat0, lon0, deg: bool = True) -> tuple:
    """
    Convert an NWU vector into ECEF coordinates.
    """

    return enu2ecefv(-w, n, u, lat0, lon0, deg=deg)


def ecef2nwu(
    x,
    y,
    z,
    lat0,
    lon0,
    h0,
    ell: Ellipsoid | None = None,
    deg: bool = True,
) -> tuple:
    """
    Convert target ECEF coordinates to North, West, Up.
    """

    e, n, u = ecef2enu(x, y, z, lat0, lon0, h0, ell, deg=deg)
    return n, -e, u


def ecef2nwuv(x, y, z, lat0, lon0, deg: bool = True) -> tuple:
    """
    Convert an ECEF vector into North, West, Up.
    """

    e, n, u = ecef2enuv(x, y, z, lat0, lon0, deg=deg)
    return n, -e, u


def geodetic2nwu(
    lat,
    lon,
    h,
    lat0,
    lon0,
    h0,
    ell: Ellipsoid | None = None,
    deg: bool = True,
) -> tuple:
    """
    Convert target geodetic coordinates to North, West, Up.
    """

    e, n, u = geodetic2enu(lat, lon, h, lat0, lon0, h0, ell, deg=deg)
    return n, -e, u


def ecef2nwu_matrix(lat0, lon0, deg: bool = True):
    """
    Rotation matrix that maps ECEF vectors into NWU coordinates.
    """

    return _ecef2nwu_rotation(lat0, lon0, deg=deg)


def nwu2ecef_matrix(lat0, lon0, deg: bool = True):
    """
    Rotation matrix that maps NWU vectors into ECEF coordinates.
    """

    return _nwu2ecef_rotation(lat0, lon0, deg=deg)
