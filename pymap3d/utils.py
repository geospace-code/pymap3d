"""all functions assume radians"""
import typing
from .ellipsoid import Ellipsoid

try:
    from numpy import hypot, cos, sin, arctan2 as atan2, radians, pi, asarray
except ImportError:
    from math import atan2, hypot, cos, sin, radians, pi

    asarray = None

__all__ = ["cart2pol", "pol2cart", "cart2sph", "sph2cart"]


def cart2pol(x: float, y: float) -> typing.Tuple[float, float]:
    """Transform Cartesian to polar coordinates"""
    return atan2(y, x), hypot(x, y)


def pol2cart(theta: float, rho: float) -> typing.Tuple[float, float]:
    """Transform polar to Cartesian coordinates"""
    return rho * cos(theta), rho * sin(theta)


def cart2sph(x: float, y: float, z: float) -> typing.Tuple[float, float, float]:
    """Transform Cartesian to spherical coordinates"""
    hxy = hypot(x, y)
    r = hypot(hxy, z)
    el = atan2(z, hxy)
    az = atan2(y, x)
    return az, el, r


def sph2cart(az: float, el: float, r: float) -> typing.Tuple[float, float, float]:
    """Transform spherical to Cartesian coordinates"""
    rcos_theta = r * cos(el)
    x = rcos_theta * cos(az)
    y = rcos_theta * sin(az)
    z = r * sin(el)
    return x, y, z


def sanitize(lat: float, ell: Ellipsoid, deg: bool) -> typing.Tuple[float, typing.Any]:
    if ell is None:
        ell = Ellipsoid()
    if asarray is not None:
        lat = asarray(lat)
    if deg:
        lat = radians(lat)

    try:
        if (abs(lat) > pi / 2).any():
            raise ValueError("-pi <= latitude <= pi")
    except (AttributeError, TypeError):
        if abs(lat) > pi / 2:
            raise ValueError("-pi <= latitude <= pi")

    return lat, ell
