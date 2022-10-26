# Copyright (c) 2014-2018 Michael Hirsch, Ph.D.
""" manipulations of sidereal time """
from datetime import datetime
from math import tau

from .timeconv import str2dt

try:
    import astropy.units as u
    from astropy.coordinates import Longitude
    from astropy.time import Time
except ImportError:
    pass


__all__ = ["datetime2sidereal", "juliandate", "greenwichsrt"]


def datetime2sidereal(time: datetime, lon_radians: float) -> float:
    """
    Convert ``datetime`` to local sidereal time

    from D. Vallado "Fundamentals of Astrodynamics and Applications"


    time : datetime.datetime
        time to convert
    lon_radians : float
        longitude (radians)

    Results
    -------

    tsr : float
        Local sidereal time
    """
    if isinstance(time, (tuple, list)):
        return [datetime2sidereal(t, lon_radians) for t in time]

    try:
        tsr = (
            Time(time)
            .sidereal_time(kind="apparent", longitude=Longitude(lon_radians, unit=u.radian))
            .radian
        )
    except NameError:
        jd = juliandate(str2dt(time))
        # %% Greenwich Sidereal time RADIANS
        gst = greenwichsrt(jd)
        # %% Algorithm 15 p. 188 rotate GST to LOCAL SIDEREAL TIME
        tsr = gst + lon_radians

    return tsr


def juliandate(time: datetime) -> float:
    """
    Python datetime to Julian time (days since Jan 1, 4713 BCE)

    from D.Vallado Fundamentals of Astrodynamics and Applications p.187
     and J. Meeus Astronomical Algorithms 1991 Eqn. 7.1 pg. 61

    Parameters
    ----------

    time : datetime.datetime
        time to convert

    Results
    -------

    jd : float
        Julian date (days since Jan 1, 4713 BCE)
    """
    if isinstance(time, (tuple, list)):
        return list(map(juliandate, time))

    if time.month < 3:
        year = time.year - 1
        month = time.month + 12
    else:
        year = time.year
        month = time.month

    A = int(year / 100.0)
    B = 2 - A + int(A / 4.0)
    C = ((time.second / 60.0 + time.minute) / 60.0 + time.hour) / 24.0

    return int(365.25 * (year + 4716)) + int(30.6001 * (month + 1)) + time.day + B - 1524.5 + C


def greenwichsrt(Jdate: float) -> float:
    """
    Convert Julian time to sidereal time

    D. Vallado Ed. 4

    Parameters
    ----------

    Jdate: float
        Julian date (since Jan 1, 4713 BCE)

    Results
    -------

    tsr : float
        Sidereal time
    """
    if isinstance(Jdate, (tuple, list)):
        return list(map(greenwichsrt, Jdate))

    # %% Vallado Eq. 3-42 p. 184, Seidelmann 3.311-1
    tUT1 = (Jdate - 2451545.0) / 36525.0

    # Eqn. 3-47 p. 188
    gmst_sec = (
        67310.54841
        + (876600 * 3600 + 8640184.812866) * tUT1
        + 0.093104 * tUT1**2
        - 6.2e-6 * tUT1**3
    )

    # 1/86400 and %(2*pi) implied by units of radians
    return gmst_sec * tau / 86400.0 % tau
