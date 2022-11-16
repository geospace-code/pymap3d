"""Utility functions

all assume radians"""

from __future__ import annotations

from math import pi
from typing import Any, Sequence, overload

try:
    from numpy import asarray
    from numpy.typing import NDArray
except ImportError:
    pass

from ._types import ArrayLike
from .ellipsoid import Ellipsoid
from .mathfun import atan2, cos, hypot, radians, sin

__all__ = ["cart2pol", "pol2cart", "cart2sph", "sph2cart"]


@overload
def cart2pol(x: float, y: float) -> tuple[float, float]:
    pass


@overload
def cart2pol(x: ArrayLike, y: ArrayLike) -> tuple[NDArray[Any], NDArray[Any]]:
    pass


def cart2pol(
    x: float | ArrayLike, y: float | ArrayLike
) -> tuple[float, float] | tuple[NDArray[Any], NDArray[Any]]:
    """Transform Cartesian to polar coordinates"""
    return atan2(y, x), hypot(x, y)


@overload
def pol2cart(theta: float, rho: float) -> tuple[float, float]:
    pass


@overload
def pol2cart(theta: ArrayLike, rho: ArrayLike) -> tuple[NDArray[Any], NDArray[Any]]:
    pass


def pol2cart(
    theta: float | ArrayLike, rho: float | ArrayLike
) -> tuple[float, float] | tuple[NDArray[Any], NDArray[Any]]:
    """Transform polar to Cartesian coordinates"""
    return rho * cos(theta), rho * sin(theta)


@overload
def cart2sph(x: float, y: float, z: float) -> tuple[float, float, float]:
    pass


@overload
def cart2sph(
    x: ArrayLike, y: ArrayLike, z: ArrayLike
) -> tuple[NDArray[Any], NDArray[Any], NDArray[Any]]:
    pass


def cart2sph(
    x: float | ArrayLike, y: float | ArrayLike, z: float | ArrayLike
) -> tuple[float, float, float] | tuple[NDArray[Any], NDArray[Any], NDArray[Any]]:
    """Transform Cartesian to spherical coordinates"""
    hxy = hypot(x, y)
    r = hypot(hxy, z)
    el = atan2(z, hxy)
    az = atan2(y, x)
    return az, el, r


@overload
def sph2cart(az: float, el: float, r: float) -> tuple[float, float, float]:
    pass


@overload
def sph2cart(
    az: ArrayLike, el: ArrayLike, r: ArrayLike
) -> tuple[NDArray[Any], NDArray[Any], NDArray[Any]]:
    pass


def sph2cart(
    az: float | ArrayLike, el: float | ArrayLike, r: float | ArrayLike
) -> tuple[float, float, float] | tuple[NDArray[Any], NDArray[Any], NDArray[Any]]:
    """Transform spherical to Cartesian coordinates"""
    rcos_theta = r * cos(el)
    x = rcos_theta * cos(az)
    y = rcos_theta * sin(az)
    z = r * sin(el)
    return x, y, z


@overload
def sanitize(lat: float, ell: Ellipsoid | None, deg: bool) -> tuple[float, Ellipsoid]:
    pass


@overload
def sanitize(lat: ArrayLike, ell: Ellipsoid | None, deg: bool) -> tuple[NDArray[Any], Ellipsoid]:
    pass


def sanitize(
    lat: float | ArrayLike, ell: Ellipsoid | None, deg: bool
) -> tuple[float | NDArray[Any], Ellipsoid]:

    if ell is None:
        ell = Ellipsoid.from_name("wgs84")

    try:
        lat = asarray(lat)
    except NameError:
        pass
    assert not isinstance(lat, Sequence)

    if deg:
        lat = radians(lat)

    try:
        if (abs(lat) > pi / 2).any():  # type: ignore[attr-defined, operator]
            raise ValueError("-pi/2 <= latitude <= pi/2")
    except AttributeError:
        if abs(lat) > pi / 2:  # type: ignore[operator]
            raise ValueError("-pi/2 <= latitude <= pi/2")

    return lat, ell
