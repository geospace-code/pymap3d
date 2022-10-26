"""
Azimuth / elevation <==> Right ascension, declination
"""

from __future__ import annotations

from datetime import datetime

from .timeconv import str2dt  # astropy can't handle xarray times (yet)
from .vallado import azel2radec as vazel2radec
from .vallado import radec2azel as vradec2azel

try:
    from astropy import units as u
    from astropy.coordinates import ICRS, AltAz, Angle, EarthLocation, SkyCoord
    from astropy.time import Time
except ImportError:
    pass

__all__ = ["radec2azel", "azel2radec"]


def azel2radec(
    az_deg: float,
    el_deg: float,
    lat_deg: float,
    lon_deg: float,
    time: datetime,
) -> tuple[float, float]:
    """
    viewing angle (az, el) to sky coordinates (ra, dec)

    Parameters
    ----------
    az_deg : float
         azimuth [degrees clockwize from North]
    el_deg : float
             elevation [degrees above horizon (neglecting aberration)]
    lat_deg : float
              observer latitude [-90, 90]
    lon_deg : float
              observer longitude [-180, 180] (degrees)
    time : datetime.datetime or str
           time of observation

    Returns
    -------
    ra_deg : float
         ecliptic right ascension (degress)
    dec_deg : float
         ecliptic declination (degrees)
    """

    try:

        obs = EarthLocation(lat=lat_deg * u.deg, lon=lon_deg * u.deg)

        direc = AltAz(
            location=obs, obstime=Time(str2dt(time)), az=az_deg * u.deg, alt=el_deg * u.deg
        )

        sky = SkyCoord(direc.transform_to(ICRS()))

        return sky.ra.deg, sky.dec.deg
    except NameError:
        return vazel2radec(az_deg, el_deg, lat_deg, lon_deg, time)


def radec2azel(
    ra_deg: float,
    dec_deg: float,
    lat_deg: float,
    lon_deg: float,
    time: datetime,
) -> tuple[float, float]:
    """
    sky coordinates (ra, dec) to viewing angle (az, el)

    Parameters
    ----------
    ra_deg : float
         ecliptic right ascension (degress)
    dec_deg : float
         ecliptic declination (degrees)
    lat_deg : float
              observer latitude [-90, 90]
    lon_deg : float
              observer longitude [-180, 180] (degrees)
    time : datetime.datetime or str
           time of observation

    Returns
    -------
    az_deg : float
             azimuth [degrees clockwize from North]
    el_deg : float
             elevation [degrees above horizon (neglecting aberration)]
    """

    try:
        obs = EarthLocation(lat=lat_deg * u.deg, lon=lon_deg * u.deg)
        points = SkyCoord(Angle(ra_deg, unit=u.deg), Angle(dec_deg, unit=u.deg), equinox="J2000.0")
        altaz = points.transform_to(AltAz(location=obs, obstime=Time(str2dt(time))))

        return altaz.az.degree, altaz.alt.degree
    except NameError:
        return vradec2azel(ra_deg, dec_deg, lat_deg, lon_deg, time)
