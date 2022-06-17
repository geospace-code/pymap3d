"""
import from Numpy, and if not available fallback to math stdlib
"""

try:
    from numpy import (
        sin,
        cos,
        sqrt,
        exp,
        log,
        inf,
        isnan,
        radians,
        tan,
        arctan as atan,
        hypot,
        degrees,
        arctan2 as atan2,
        arcsin as asin,
        arcsinh as asinh,
        arctanh as atanh,
        power,
    )
except ImportError:
    from math import sin, cos, sqrt, exp, log, inf, isnan, radians, tan, atan, hypot, degrees, atan2, asin, asinh, atanh  # type: ignore

    def power(x, y):  # type: ignore
        return pow(x, y)
