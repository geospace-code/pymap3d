"""
import from Numpy, and if not available fallback to math stdlib
"""

from __future__ import annotations

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
        isclose,
        isnan,
        linspace,
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
        isclose,
        isnan,
        log,
        radians,
        sin,
        sqrt,
        tan,
    )

    def linspace(start: float, stop: float, num: int) -> list[float]:  # type: ignore
        """
        create a list of "num" evenly spaced numbers using range and increment,
        including endpoint "stop"
        """
        step = (stop - start) / (num - 1)
        return [start + i * step for i in range(num)]

    def power(x, y):  # type: ignore
        return pow(x, y)

    def sign(x) -> float:  # type: ignore
        """signum"""
        if x < 0:
            y = -1.0
        elif x > 0:
            y = 1.0
        else:
            y = 0.0

        return y

    try:
        import math.cbrt as cbrt  # type: ignore
    except ImportError:

        def cbrt(x) -> float:  # type: ignore
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
    "isclose",
    "isnan",
    "log",
    "power",
    "radians",
    "sign",
    "sin",
    "sqrt",
    "tan",
]
