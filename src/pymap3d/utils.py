"""Utility functions

all assume radians"""

from __future__ import annotations
import typing

from .ellipsoid import Ellipsoid

try:
    from numpy import hypot, cos, sin, arctan2 as atan2, radians, pi, asarray, ndarray
except ImportError:
    from math import atan2, hypot, cos, sin, radians, pi  # type: ignore

    asarray = None
    ndarray = typing.Any  # type: ignore

__all__ = ["cart2pol", "pol2cart", "cart2sph", "sph2cart", "sign"]


def sign(x: float) -> float:
    """ signum function """
    if x < 0:
        y = -1.0
    elif x > 0:
        y = 1.0
    else:
        y = 0.0

    return y


def cart2pol(x: float, y: float) -> tuple[float, float]:
    """Transform Cartesian to polar coordinates"""
    return atan2(y, x), hypot(x, y)


def pol2cart(theta: float, rho: float) -> tuple[float, float]:
    """Transform polar to Cartesian coordinates"""
    return rho * cos(theta), rho * sin(theta)


def cart2sph(
    x: float | ndarray, y: float | ndarray, z: float | ndarray
) -> tuple[float | ndarray, float | ndarray, float | ndarray]:
    """Transform Cartesian to spherical coordinates"""
    hxy = hypot(x, y)
    r = hypot(hxy, z)
    el = atan2(z, hxy)
    az = atan2(y, x)
    return az, el, r


def sph2cart(
    az: float | ndarray, el: float | ndarray, r: float | ndarray
) -> tuple[float | ndarray, float | ndarray, float | ndarray]:
    """Transform spherical to Cartesian coordinates"""
    rcos_theta = r * cos(el)
    x = rcos_theta * cos(az)
    y = rcos_theta * sin(az)
    z = r * sin(el)
    return x, y, z


def sanitize(lat: float | ndarray, ell: Ellipsoid, deg: bool) -> tuple[float | ndarray, Ellipsoid]:
    if ell is None:
        ell = Ellipsoid()
    if asarray is not None:
        lat = asarray(lat)
    if deg:
        lat = radians(lat)

    if asarray is not None:
        if (abs(lat) > pi / 2).any():
            raise ValueError("-pi/2 <= latitude <= pi/2")
    else:
        if abs(lat) > pi / 2:
            raise ValueError("-pi/2 <= latitude <= pi/2")

    return lat, ell
