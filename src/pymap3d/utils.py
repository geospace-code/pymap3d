"""Utility functions

all assume radians"""

from __future__ import annotations

from ._typing import FloatLike, FloatArray

from .mathfun import atan2, cos, hypot, sin

__all__ = ["cart2pol", "pol2cart", "cart2sph", "sph2cart"]


def cart2pol(x: FloatLike, y: FloatLike) -> tuple:
    """Transform Cartesian to polar coordinates

    Parameters
    ----------

    x : array-like float
        x coordinate
    y : array-like float
        y coordinate

    Returns
    -------

    theta : array-like float
        angle in radians
    rho : array-like float
        radius
    """
    return atan2(y, x), hypot(x, y)


def pol2cart(theta: FloatLike, rho: FloatLike) -> tuple:
    """Transform polar to Cartesian coordinates

    Parameters
    ----------

    theta : array-like float
        angle in radians
    rho : array-like float
        radius

    Returns
    -------

    x : array-like float
        x coordinate
    y : array-like float
        y coordinate
    """
    return rho * cos(theta), rho * sin(theta)


def cart2sph(x: FloatLike, y: FloatLike, z: FloatLike) -> tuple:
    """Transform Cartesian to spherical coordinates

    Parameters
    ----------

    x : array-like float
        x coordinate
    y : array-like float
        y coordinate
    z : array-like float
        z coordinate

    Returns
    -------

    az : array-like float
        azimuth angle in radians
    el : array-like float
        elevation angle in radians
    r : array-like float
        radius
    """
    hxy = hypot(x, y)
    r = hypot(hxy, z)
    el = atan2(z, hxy)
    az = atan2(y, x)
    return az, el, r


def sph2cart(az: FloatArray, el: FloatArray, r: FloatLike) -> tuple:
    """Transform spherical to Cartesian coordinates

    Parameters
    ----------

    az : array-like float
        azimuth angle in radians
    el : array-like float
        elevation angle in radians
    r : array-like float
        radius

    Returns
    -------

    x : array-like float
        x coordinate
    y : array-like float
        y coordinate
    z : array-like float
        z coordinate
    """
    rcos_theta = r * cos(el)
    x = rcos_theta * cos(az)
    y = rcos_theta * sin(az)
    z = r * sin(el)
    return x, y, z
