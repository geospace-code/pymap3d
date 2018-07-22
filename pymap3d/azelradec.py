"""
Azimuth / elevation <==> Right ascension, declination
"""
from typing import Tuple, Union
from datetime import datetime
import numpy as np
from .timeconv import str2dt
try:
    from astropy.time import Time
    from astropy import units as u
    from astropy.coordinates import Angle, SkyCoord, EarthLocation, AltAz, ICRS
except ImportError:
    from .vallado import vazel2radec, vradec2azel
    Time = None


def azel2radec(az_deg: float, el_deg: float,
               lat_deg: float, lon_deg: float,
               time: Union[str, datetime]) -> Tuple[float, float]:
    """convert astronomical target horizontal azimuth, elevation to
       ecliptic right ascension, declination (degrees)
    """

    if Time is None:  # non-AstroPy method, less accurate
        return vazel2radec(az_deg, el_deg, lat_deg, lon_deg, time)

    t = str2dt(time)

    obs = EarthLocation(lat=lat_deg * u.deg, lon=lon_deg * u.deg)

    direc = AltAz(location=obs, obstime=Time(t),
                  az=az_deg * u.deg, alt=el_deg * u.deg)

    sky = SkyCoord(direc.transform_to(ICRS()))

    return sky.ra.deg, sky.dec.deg


def radec2azel(ra_deg: float, dec_deg: float,
               lat_deg: float, lon_deg: float,
               time: Union[str, datetime]) -> Tuple[float, float]:
    """convert astronomical target ecliptic right ascension, declination to
       horizontal azimuth, eelvation (degrees)
    """
    if Time is None:
        return vradec2azel(ra_deg, dec_deg, lat_deg, lon_deg, time)
# %% input trapping
    t = str2dt(time)
    lat = np.atleast_1d(lat_deg)
    lon = np.atleast_1d(lon_deg)
    ra = np.atleast_1d(ra_deg)
    dec = np.atleast_1d(dec_deg)

    if not(lat.size == 1 & lon.size == 1):
        raise ValueError('radec2azel is designed for one observer and one or more points (ra,dec).')

    if ra.shape != dec.shape:
        raise ValueError('ra and dec must be the same shape ndarray')

    obs = EarthLocation(lat=lat * u.deg,
                        lon=lon * u.deg)

    points = SkyCoord(Angle(ra, unit=u.deg),
                      Angle(dec, unit=u.deg),
                      equinox='J2000.0')

    altaz = points.transform_to(AltAz(location=obs, obstime=Time(t)))

    return altaz.az.degree, altaz.alt.degree
