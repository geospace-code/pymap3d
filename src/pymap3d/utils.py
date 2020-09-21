"""Utility functions

all assume radians"""
import typing
from .ellipsoid import Ellipsoid

try:
    from numpy import hypot, cos, sin, arctan2 as atan2, radians, pi, asarray
except ImportError:
    from math import atan2, hypot, cos, sin, radians, pi

    asarray = None

try:
    from numpy.typing import ArrayLike
except ImportError:
    ArrayLike = typing.Any

__all__ = ["cart2pol", "pol2cart", "cart2sph", "sph2cart", "sign"]


def sign(x: ArrayLike) -> ArrayLike:
    """ signum function """
    if x < 0:
        y = -1.0
    elif x > 0:
        y = 1.0
    else:
        y = 0.0

    return y


def cart2pol(x: ArrayLike, y: ArrayLike) -> typing.Tuple[ArrayLike, ArrayLike]:
    """Transform Cartesian to polar coordinates"""
    return atan2(y, x), hypot(x, y)


def pol2cart(theta: ArrayLike, rho: ArrayLike) -> typing.Tuple[ArrayLike, ArrayLike]:
    """Transform polar to Cartesian coordinates"""
    return rho * cos(theta), rho * sin(theta)


def cart2sph(x: ArrayLike, y: ArrayLike, z: ArrayLike) -> typing.Tuple[ArrayLike, ArrayLike, ArrayLike]:
    """Transform Cartesian to spherical coordinates"""
    hxy = hypot(x, y)
    r = hypot(hxy, z)
    el = atan2(z, hxy)
    az = atan2(y, x)
    return az, el, r


def sph2cart(az: ArrayLike, el: ArrayLike, r: ArrayLike) -> typing.Tuple[ArrayLike, ArrayLike, ArrayLike]:
    """Transform spherical to Cartesian coordinates"""
    rcos_theta = r * cos(el)
    x = rcos_theta * cos(az)
    y = rcos_theta * sin(az)
    z = r * sin(el)
    return x, y, z


def sanitize(lat: ArrayLike, ell: Ellipsoid, deg: bool) -> typing.Tuple[ArrayLike, Ellipsoid]:
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
