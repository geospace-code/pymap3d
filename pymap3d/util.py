from numpy import arctan2, hypot, cos, sin, pi
from typing import Tuple

__all__ = ['cart2pol', 'pol2cart', 'cart2sph', 'sph2cart', 'wrapToPi', 'wrapTo2Pi']

def cart2pol(x: float, y: float) -> Tuple[float, float]:
    '''Transform Cartesian to polar coordinates'''
    return arctan2(y, x), hypot(x, y)

def pol2cart(theta: float, rho: float) -> Tuple[float, float]:
    '''Transform polar to Cartesian coordinates'''
    return  rho * cos(theta), rho * sin(theta)

def cart2sph(x: float, y: float, z: float) -> Tuple[float, float, float]:
    '''Transform Cartesian to spherical coordinates'''
    hxy = hypot(x, y)
    r = hypot(hxy, z)
    el = arctan2(z, hxy)
    az = arctan2(y, x)
    return az, el, r

def sph2cart(az: float, el: float, r: float) -> Tuple[float, float, float]:
    '''Transform spherical to Cartesian coordinates'''
    rcos_theta = r * cos(el)
    x = rcos_theta * cos(az)
    y = rcos_theta * sin(az)
    z = r * sin(el)
 return x, y, z

def wrapToPi(x: float) -> float:
    '''wraps angles to [-pi, pi)'''
    return ((x + pi) % (2 * pi)) - pi

def wrapTo2Pi(x: float) -> float:
    '''wraps angles to [0, 2*pi]'''
    return x % (2 * pi)





