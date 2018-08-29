"""
Azimuth / elevation <==> Right ascension, declination
"""
from typing import Tuple
from datetime import datetime
import numpy as np
from .vallado import azel2radec as vazel2radec, radec2azel as vradec2azel
try:
    from astropy.time import Time
    from astropy import units as u
    from astropy.coordinates import Angle, SkyCoord, EarthLocation, AltAz, ICRS
except ImportError:
    Time = None


def azel2radec(az_deg: float, el_deg: float,
               lat_deg: float, lon_deg: float,
               time: datetime, usevallado: bool=False) -> Tuple[float, float]:
    """
    viewing angle (az, el) to sky coordinates (ra, dec)

    inputs
    ------
    azimuth: degrees clockwize from North
    elevation: degrees above horizon (neglecting aberration)
    observer latitude [-90, 90], longitude [-180, 180] (degrees)
    time: datetime of observation

    Outputs
    -------
    ecliptic right ascension, declination (degrees)
    """

    if usevallado or Time is None:  # non-AstroPy method, less accurate
        return vazel2radec(az_deg, el_deg, lat_deg, lon_deg, time)

    obs = EarthLocation(lat=lat_deg * u.deg, lon=lon_deg * u.deg)

    direc = AltAz(location=obs, obstime=Time(time),
                  az=az_deg * u.deg, alt=el_deg * u.deg)

    sky = SkyCoord(direc.transform_to(ICRS()))

    return sky.ra.deg, sky.dec.deg


def radec2azel(ra_deg: float, dec_deg: float,
               lat_deg: float, lon_deg: float,
               time: datetime, usevallado: bool=False) -> Tuple[float, float]:
    """
    sky coordinates (ra, dec) to viewing angle (az, el)

    inputs
    ------
    ecliptic right ascension, declination (degrees)
    observer latitude [-90, 90], longitude [-180, 180] (degrees)
    time: datetime of observation

    Outputs
    -------
    azimuth: degrees clockwize from North
    elevation: degrees above horizon (neglecting aberration)
    """

    if usevallado or Time is None:
        return vradec2azel(ra_deg, dec_deg, lat_deg, lon_deg, time)
# %% input trapping
    lat = np.atleast_1d(lat_deg)
    lon = np.atleast_1d(lon_deg)
    ra = np.atleast_1d(ra_deg)
    dec = np.atleast_1d(dec_deg)

    obs = EarthLocation(lat=lat * u.deg,
                        lon=lon * u.deg)

    points = SkyCoord(Angle(ra, unit=u.deg),
                      Angle(dec, unit=u.deg),
                      equinox='J2000.0')

    altaz = points.transform_to(AltAz(location=obs, obstime=Time(time)))

    return altaz.az.degree, altaz.alt.degree
