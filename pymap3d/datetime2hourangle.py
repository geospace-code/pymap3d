# Copyright (c) 2014-2018 Michael Hirsch, Ph.D.
from math import pi
from typing import Union, List
from datetime import datetime
from .timeconv import str2dt
#
nan = float('nan')
try:
    from astropy.time import Time
    import astropy.units as u
    from astropy.coordinates import Longitude
except ImportError:
    Time = None

"""
The "usevallado" datetime to julian runs 4 times faster than astropy.
However, AstroPy is more accurate.
"""


def datetime2sidereal(t: Union[str, datetime],
                      lon_radians: float,
                      usevallado: bool=True) -> Union[float, List[float]]:
    """
    Convert ``datetime`` to sidereal time

    :algorithm: D. Vallado Fundamentals of Astrodynamics and Applications


    t
        Python datetime

    lon
        longitude in RADIANS
    """
    if usevallado:
        jd = datetime2julian(t)
# %% Greenwich Sidereal time RADIANS
        gst = julian2sidereal(jd)
# %% Algorithm 15 p. 188 rotate GST to LOCAL SIDEREAL TIME
        tsr: Union[float, List[float]]
        if isinstance(gst, float):
            tsr = gst + lon_radians  # radians
        else:
            tsr = [g + lon_radians for g in gst]
    else:  # astropy
        if Time is not None:
            tsr = Time(t).sidereal_time(kind='apparent',
                                        longitude=Longitude(lon_radians, unit=u.radian)).radian
        else:
            raise ImportError('AstroPy required, or use "usevallado=True"')

    return tsr


def datetime2julian(time: Union[str, datetime]) -> Union[float, List[float]]:
    """
    Python datetime to Julian time

    from D.Vallado Fundamentals of Astrodynamics and Applications p.187
     and J. Meeus Astronomical Algorithms 1991 Eqn. 7.1 pg. 61
    """

    def _dt2julian(t: datetime) -> float:
        if t.month < 3:
            year = t.year - 1
            month = t.month + 12
        else:
            year = t.year
            month = t.month

        A = int(year / 100.0)
        B = 2 - A + int(A / 4.)
        C = ((t.second / 60. + t.minute) / 60. + t.hour) / 24.

        jd = (int(365.25 * (year + 4716)) +
              int(30.6001 * (month + 1)) + t.day + B - 1524.5 + C)

        return jd

    T = str2dt(time)

    if isinstance(T, datetime):
        return _dt2julian(T)
    else:
        return [_dt2julian(t) for t in T]


def julian2sidereal(juliandate: Union[float, List[float]]) -> Union[float, List[float]]:
    """
    Convert Julian time to sidereal time

    D. Vallado Ed. 4

    input:

    juliandate
        Julian centuries from J2000.0

    """

    def _jd2sr(jd: float):
        # %% Vallado Eq. 3-42 p. 184, Seidelmann 3.311-1
        tUT1 = (jd - 2451545.0) / 36525.

        # Eqn. 3-47 p. 188
        gmst_sec = (67310.54841 + (876600 * 3600 + 8640184.812866) *
                    tUT1 + 0.093104 * tUT1**2 - 6.2e-6 * tUT1**3)

        # 1/86400 and %(2*pi) implied by units of radians
        tsr = gmst_sec * (2 * pi) / 86400. % (2 * pi)

        return tsr

    if isinstance(juliandate, float):
        return _jd2sr(juliandate)
    else:
        return [_jd2sr(jd) for jd in juliandate]
