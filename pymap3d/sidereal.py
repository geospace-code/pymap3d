# Copyright (c) 2014-2018 Michael Hirsch, Ph.D.
""" manipulations of sidereal time """
from math import pi
from datetime import datetime
import numpy as np
from .timeconv import str2dt
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

__all__ = ['datetime2sidereal', 'juliandate', 'julian2sidereal']


def datetime2sidereal(time: datetime,
                      lon_radians: float,
                      usevallado: bool = True) -> float:
    """
    Convert ``datetime`` to sidereal time

    from D. Vallado "Fundamentals of Astrodynamics and Applications"


    time : datetime.datetime
        time to convert
    lon_radians : float
        longitude (radians)
    usevallado : bool, optional
        use vallado instead of AstroPy (default is Vallado)

    Results
    -------

    tsr : float
        Sidereal time
    """
    usevallado = usevallado or Time is None
    if usevallado:
        jd = juliandate(str2dt(time))
# %% Greenwich Sidereal time RADIANS
        gst = julian2sidereal(jd)
# %% Algorithm 15 p. 188 rotate GST to LOCAL SIDEREAL TIME
        tsr = gst + lon_radians
    else:
        tsr = Time(time).sidereal_time(kind='apparent',
                                       longitude=Longitude(lon_radians, unit=u.radian)).radian

    return tsr


def juliandate(time: datetime) -> float:
    """
    Python datetime to Julian time

    from D.Vallado Fundamentals of Astrodynamics and Applications p.187
     and J. Meeus Astronomical Algorithms 1991 Eqn. 7.1 pg. 61

    Parameters
    ----------

    time : datetime.datetime
        time to convert

    Results
    -------

    jd : float
        Julian date
    """

    times = np.atleast_1d(time)
    assert times.ndim == 1

    jd = np.empty(times.size)
    for i, t in enumerate(times):
        if t.month < 3:
            year = t.year - 1
            month = t.month + 12
        else:
            year = t.year
            month = t.month

        A = int(year / 100.0)
        B = 2 - A + int(A / 4.)
        C = ((t.second / 60. + t.minute) / 60. + t.hour) / 24.

        jd[i] = (int(365.25 * (year + 4716)) +
                 int(30.6001 * (month + 1)) + t.day + B - 1524.5 + C)

    return jd.squeeze()


def julian2sidereal(Jdate: float) -> float:
    """
    Convert Julian time to sidereal time

    D. Vallado Ed. 4

    Parameters
    ----------

    Jdate: float
        Julian centuries from J2000.0

    Results
    -------

    tsr : float
        Sidereal time

    """

    jdate = np.atleast_1d(Jdate)
    assert jdate.ndim == 1

    tsr = np.empty(jdate.size)
    for i, jd in enumerate(jdate):
        # %% Vallado Eq. 3-42 p. 184, Seidelmann 3.311-1
        tUT1 = (jd - 2451545.0) / 36525.

        # Eqn. 3-47 p. 188
        gmst_sec = (67310.54841 + (876600 * 3600 + 8640184.812866) *
                    tUT1 + 0.093104 * tUT1**2 - 6.2e-6 * tUT1**3)

        # 1/86400 and %(2*pi) implied by units of radians
        tsr[i] = gmst_sec * (2 * pi) / 86400. % (2 * pi)

    return tsr.squeeze()
