"""
import from Numpy, and if not available fallback to math stdlib
"""

try:
    from numpy import arcsin as asin
    from numpy import arcsinh as asinh
    from numpy import arctan as atan
    from numpy import arctan2 as atan2
    from numpy import arctanh as atanh
    from numpy import (
        cbrt,
        cos,
        degrees,
        exp,
        hypot,
        inf,
        isnan,
        log,
        power,
        radians,
        sign,
        sin,
        sqrt,
        tan,
    )
except ImportError:
    from math import (  # type: ignore
        asin,
        asinh,
        atan,
        atan2,
        atanh,
        cos,
        degrees,
        exp,
        hypot,
        inf,
        isnan,
        log,
        radians,
        sin,
        sqrt,
        tan,
    )

    def power(x, y):  # type: ignore
        return pow(x, y)

    def sign(x) -> float:  # type: ignore
        """signum function"""
        if x < 0:
            y = -1.0
        elif x > 0:
            y = 1.0
        else:
            y = 0.0

        return y

    def cbrt(x) -> float:  # type: ignore
        """math.cbrt was added in Python 3.11"""
        return x ** (1 / 3)


__all__ = [
    "asin",
    "asinh",
    "atan",
    "atan2",
    "atanh",
    "cbrt",
    "cos",
    "degrees",
    "exp",
    "hypot",
    "inf",
    "isnan",
    "log",
    "power",
    "radians",
    "sign",
    "sin",
    "sqrt",
    "tan",
]
