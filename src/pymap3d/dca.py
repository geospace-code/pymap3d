"""
Transforms involving DCA (Downrange, Crossrange, Above)

This module provides functions to convert coordinates between DCA and
other coordinate systems such as ENU (East, North, Up), NED (North, East, Down),
ECEF (Earth-Centered, Earth-Fixed), and geodetic coordinates (latitude, longitude, altitude).

It also includes transformations to/from AER (Azimuth, Elevation, Range).

Functions:
----------
- enu2dca: Convert ENU to DCA coordinates.
- dca2enu: Convert DCA to ENU coordinates.
- dca2ned: Convert DCA to NED coordinates.
- ned2dca: Convert NED to DCA coordinates.
- ecef2dca: Convert ECEF to DCA coordinates.
- dca2ecef: Convert DCA to ECEF coordinates.
- geodetic2dca: Convert geodetic to DCA coordinates.
- dca2geodetic: Convert DCA to geodetic coordinates.
- aer2dca: Convert AER to DCA coordinates.
- dca2aer: Convert DCA to AER coordinates.
"""

from __future__ import annotations

from .mathfun import sin, cos, radians
from .ecef import ecef2enu, enu2ecef
from .enu import geodetic2enu, enu2geodetic, aer2enu, enu2aer
from .ellipsoid import Ellipsoid

__all__ = [
    "enu2dca",
    "dca2enu",
    "dca2ned",
    "ned2dca",
    "ecef2dca",
    "dca2ecef",
    "geodetic2dca",
    "dca2geodetic",
    "aer2dca",
    "dca2aer",
]


def enu2dca(e, n, u, heading, deg: bool = True):
    """
    Converts ENU (East, North, Up) coordinates to DCA (Downrange, Crossrange, Above).
    """

    if deg:
        heading = radians(heading)

    dr = e * sin(heading) + n * cos(heading)
    cr = -e * cos(heading) + n * sin(heading)

    return dr, cr, u


def dca2enu(dr, cr, above, heading, deg: bool = True):
    """
    Converts DCA (Downrange, Crossrange, Above) coordinates to ENU (East, North, Up).
    """

    if deg:
        heading = radians(heading)

    e = dr * sin(heading) - cr * cos(heading)
    n = dr * cos(heading) + cr * sin(heading)

    return e, n, above


def dca2ned(dr, cr, above, heading, deg: bool = True):
    """
    Converts DCA (Downrange, Crossrange, Above) coordinates to NED (North, East, Down).
    """
    e, n, u = dca2enu(dr, cr, above, heading, deg=deg)
    return n, e, -u


def ned2dca(n, e, d, heading, deg: bool = True):
    """
    Converts NED (North, East, Down) coordinates to DCA (Downrange, Crossrange, Above).
    """
    dr, cr, above = enu2dca(e, n, -d, heading, deg=deg)
    return dr, cr, above


def ecef2dca(x, y, z, lat0, lon0, h0, heading, ell: Ellipsoid | None = None, deg: bool = True):
    """
    Converts ECEF (Earth-Centered, Earth-Fixed) coordinates to DCA (Downrange, Crossrange, Above).
    """

    e, n, u = ecef2enu(x, y, z, lat0, lon0, h0, ell, deg=deg)
    return enu2dca(e, n, u, heading, deg=deg)


def dca2ecef(
    dr,
    cr,
    above,
    lat0,
    lon0,
    h0,
    heading,
    ell: Ellipsoid | None = None,
    deg: bool = True,
):
    """
    Converts DCA (Downrange, Crossrange, Above) coordinates to ECEF (Earth-Centered, Earth-Fixed) coordinates.
    """

    e, n, u = dca2enu(dr, cr, above, heading, deg=deg)
    return enu2ecef(e, n, u, lat0, lon0, h0, ell, deg=deg)


def geodetic2dca(
    lat, lon, h, lat0, lon0, h0, heading, ell: Ellipsoid | None = None, deg: bool = True
):
    """
    Converts geodetic coordinates (latitude, longitude, altitude) to DCA (Downrange, Crossrange, Above) coordinates.
    """
    e, n, u = geodetic2enu(lat, lon, h, lat0, lon0, h0, ell, deg=deg)
    return enu2dca(e, n, u, heading, deg=deg)


def dca2geodetic(
    dr,
    cr,
    above,
    lat0,
    lon0,
    h0,
    heading,
    ell: Ellipsoid | None = None,
    deg: bool = True,
):
    """
    Converts DCA (Downrange, Crossrange, Above) coordinates to geodetic coordinates (latitude, longitude, altitude).
    """
    e, n, u = dca2enu(dr, cr, above, heading, deg=deg)
    return enu2geodetic(e, n, u, lat0, lon0, h0, ell, deg=deg)


def aer2dca(az, el, srange, heading, deg: bool = True):
    """
    Converts AER (Azimuth, Elevation, Range) coordinates to DCA (Downrange, Crossrange, Above).
    """
    e, n, u = aer2enu(az, el, srange, deg=deg)
    return enu2dca(e, n, u, heading, deg=deg)


def dca2aer(dr, cr, above, heading, deg: bool = True):
    """
    Converts DCA (Downrange, Crossrange, Above) coordinates to AER (Azimuth, Elevation, Range).
    """
    e, n, u = dca2enu(dr, cr, above, heading, deg=deg)
    return enu2aer(e, n, u, deg=deg)
