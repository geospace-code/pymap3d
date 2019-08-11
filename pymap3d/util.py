from numpy import arctan2, hypot, cos, sin, pi
import numpy as np
from typing import Tuple

__all__ = ['cart2pol', 'pol2cart', 'cart2sph', 'sph2cart', 'wrapToPi', 'wrapTo2Pi']

def cart2pol(x: float, y: float) -> Tuple[float, float]:
    '''Transform Cartesian to polar coordinates'''
    return np.arctan2(y, x), np.hypot(x, y)

def pol2cart(theta: float, rho: float) -> Tuple[float, float]:
    '''Transform polar to Cartesian coordinates'''
    return  rho * np.cos(theta), rho * np.sin(theta)

def cart2sph(x: float, y: float, z: float) -> Tuple[float, float, float]:
    '''Transform Cartesian to spherical coordinates'''
    hxy = np.hypot(x, y)
    r = np.hypot(hxy, z)
    el = np.arctan2(z, hxy)
    az = np.arctan2(y, x)
    return az, el, r

def sph2cart(az: float, el: float, r: float) -> Tuple[float, float, float]:
    '''Transform spherical to Cartesian coordinates'''
    rcos_theta = r * np.cos(el)
    x = rcos_theta * np.cos(az)
    y = rcos_theta * np.sin(az)
    z = r * np.sin(el)
    return x, y, z

def wrapToPi(x: float) -> float:
    '''wraps angles to [-pi, pi)'''
    return ((x + pi) % (2 * pi)) - pi

def wrapTo2Pi(x: float) -> float:
    '''wraps angles to [0, 2*pi]'''
    return x % (2 * pi)





