"""Transforms involving SEZ South East Zenith."""

from __future__ import annotations

from .ecef import ecef2enu, ecef2enuv, enu2ecef
from .ellipsoid import Ellipsoid
from .enu import aer2enu, enu2aer, enu2ecefv, geodetic2enu, enu2geodetic
from .frames import _ecef2sez_rotation, _sez2ecef_rotation

__all__ = [
    "aer2sez",
    "sez2aer",
    "sez2geodetic",
    "sez2ecef",
    "sez2ecefv",
    "ecef2sez",
    "ecef2sezv",
    "geodetic2sez",
    "ecef2sez_matrix",
    "sez2ecef_matrix",
]


def aer2sez(az, elev, slantRange, deg: bool = True) -> tuple:
    """
    Convert azimuth, elevation, range to South, East, Zenith.
    """

    e, n, u = aer2enu(az, elev, slantRange, deg=deg)
    return -n, e, u


def sez2aer(s, e, z, deg: bool = True) -> tuple:
    """
    Convert South, East, Zenith to azimuth, elevation, range.
    """

    return enu2aer(e, -s, z, deg=deg)


def sez2geodetic(
    s,
    e,
    z,
    lat0,
    lon0,
    h0,
    ell: Ellipsoid | None = None,
    deg: bool = True,
) -> tuple:
    """
    Convert South, East, Zenith to target geodetic coordinates.
    """

    return enu2geodetic(e, -s, z, lat0, lon0, h0, ell, deg=deg)


def sez2ecef(
    s,
    e,
    z,
    lat0,
    lon0,
    h0,
    ell: Ellipsoid | None = None,
    deg: bool = True,
) -> tuple:
    """
    Convert South, East, Zenith to target ECEF coordinates.
    """

    return enu2ecef(e, -s, z, lat0, lon0, h0, ell, deg=deg)


def sez2ecefv(s, e, z, lat0, lon0, deg: bool = True) -> tuple:
    """
    Convert a SEZ vector into ECEF coordinates.
    """

    return enu2ecefv(e, -s, z, lat0, lon0, deg=deg)


def ecef2sez(
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
    Convert target ECEF coordinates to South, East, Zenith.
    """

    e, n, u = ecef2enu(x, y, z, lat0, lon0, h0, ell, deg=deg)
    return -n, e, u


def ecef2sezv(x, y, z, lat0, lon0, deg: bool = True) -> tuple:
    """
    Convert an ECEF vector into South, East, Zenith.
    """

    e, n, u = ecef2enuv(x, y, z, lat0, lon0, deg=deg)
    return -n, e, u


def geodetic2sez(
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
    Convert target geodetic coordinates to South, East, Zenith.
    """

    e, n, u = geodetic2enu(lat, lon, h, lat0, lon0, h0, ell, deg=deg)
    return -n, e, u


def ecef2sez_matrix(lat0, lon0, deg: bool = True):
    """
    Rotation matrix that maps ECEF vectors into SEZ coordinates.
    """

    return _ecef2sez_rotation(lat0, lon0, deg=deg)


def sez2ecef_matrix(lat0, lon0, deg: bool = True):
    """
    Rotation matrix that maps SEZ vectors into ECEF coordinates.
    """

    return _sez2ecef_rotation(lat0, lon0, deg=deg)
