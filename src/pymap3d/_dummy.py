from datetime import datetime


def eci2ecef(x,
    y,
    z,
    time: datetime,
    force_non_astropy: bool = False,
    *,
    delta_ut1: float = 0.0,
    xp: float = 0.0,
    yp: float = 0.0):
    raise ImportError("Numpy required for eci2ecef")


def ecef2eci(x,
    y,
    z,
    time: datetime,
    force_non_astropy: bool = False,
    *,
    delta_ut1: float = 0.0,
    xp: float = 0.0,
    yp: float = 0.0,):
    raise ImportError("Numpy required for ecef2eci")
