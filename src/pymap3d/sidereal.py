"""manipulations of sidereal time"""

from datetime import datetime
from math import tau
import sys
import logging

from .timeconv import str2dt
from .earth_orientation import (
    greenwich_apparent_sidereal_time,
    greenwich_mean_sidereal_time,
    juliandate as _juliandate,
    utc_to_tt,
    utc_with_offset,
)

try:
    import astropy.units as u
    from astropy.coordinates import Longitude
    from astropy.time import Time
except ImportError:
    pass


__all__ = ["datetime2sidereal", "juliandate", "greenwichsrt"]


def datetime2sidereal(
    time: datetime,
    lon_radians: float,
    force_non_astropy: bool = False,
    *,
    delta_ut1: float = 0.0,
):
    """
    Convert ``datetime`` to local sidereal time

    from D. Vallado "Fundamentals of Astrodynamics and Applications"


    time : datetime.datetime
        time to convert
    lon_radians : float
        longitude (radians)
    force_non_astropy : bool
        if True, force use of less accurate Numpy implementation even if Astropy is available
    delta_ut1 : float, optional
        UT1-UTC in seconds for the pure-Python path. Defaults to ``0.0``.

    Results
    -------

    tsr : float
        Local sidereal time
    """

    if isinstance(time, (tuple, list)):
        return [
            datetime2sidereal(
                t,
                lon_radians,
                force_non_astropy=force_non_astropy,
                delta_ut1=delta_ut1,
            )
            for t in time
        ]

    if "astropy" in sys.modules and not force_non_astropy:
        return datetime2sidereal_astropy(time, lon_radians)
    else:
        logging.debug(f"{__name__}: pure-Python apparent sidereal implementation")
        return datetime2sidereal_vallado(time, lon_radians, delta_ut1=delta_ut1)


def datetime2sidereal_astropy(t: datetime, lon_radians: float):
    """datetime to sidereal time using astropy
    see datetime2sidereal() for description
    """

    at = Time(t)
    tsr = at.sidereal_time(
        kind="apparent", longitude=Longitude(lon_radians, unit=u.radian)
    )
    return tsr.radian


def datetime2sidereal_vallado(
    t: datetime, lon_radians: float, *, delta_ut1: float = 0.0
):
    """datetime to sidereal time using Vallado methods
    see datetime2sidereal() for description
    """

    jd_ut1 = juliandate(utc_with_offset(str2dt(t), delta_ut1))
    jd_tt = juliandate(utc_to_tt(str2dt(t)))
    gst = greenwich_apparent_sidereal_time(jd_ut1, jd_tt)
    return (gst + lon_radians) % tau


def juliandate(time: datetime):
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

    return _juliandate(time)


def greenwichsrt(Jdate: float):
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

    return greenwich_mean_sidereal_time(Jdate)
