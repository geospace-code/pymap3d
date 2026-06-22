"""
Azimuth / elevation <==> Right ascension, declination
"""

from __future__ import annotations

from datetime import datetime
import sys

from .vallado import azel2radec as vazel2radec
from .vallado import radec2azel as vradec2azel
from ._typing import FloatLike

try:
    from astropy import units as u
    from astropy.coordinates import ICRS, AltAz, Angle, EarthLocation, SkyCoord
    from astropy.time import Time
except ImportError:
    pass

__all__ = ["radec2azel", "azel2radec"]


def azel2radec(
    az_deg: FloatLike,
    el_deg: FloatLike,
    lat_deg: FloatLike,
    lon_deg: FloatLike,
    time: datetime,
    force_non_astropy: bool = False,
    *,
    delta_ut1: float = 0.0,
) -> tuple:
    """
    viewing angle (az, el) to sky coordinates (ra, dec)

    Parameters
    ----------
    az_deg : array-like float
         azimuth [degrees clockwize from North]
    el_deg : array-like float
             elevation [degrees above horizon (neglecting aberration)]
    lat_deg : array-like float
              observer latitude [-90, 90]
    lon_deg : array-like float
              observer longitude [-180, 180] (degrees)
    time : datetime.datetime or str
           time of observation
    force_non_astropy : bool
        if True, force use of less accurate Numpy implementation even if Astropy is available

    Returns
    -------
    ra_deg : array-like float
         ecliptic right ascension (degress)
    dec_deg : array-like float
         ecliptic declination (degrees)
    """

    if force_non_astropy or "astropy" not in sys.modules:
        return vazel2radec(
            az_deg, el_deg, lat_deg, lon_deg, time, delta_ut1=delta_ut1
        )
    else:
        return azel2radec_astropy(az_deg, el_deg, lat_deg, lon_deg, time)


def azel2radec_astropy(
    az_deg: FloatLike, el_deg: FloatLike, lat_deg: FloatLike, lon_deg: FloatLike, time: datetime
) -> tuple:
    """azel2radec using Astropy
    see azel2radec() for description
    """
    obs = EarthLocation(lat=lat_deg * u.deg, lon=lon_deg * u.deg)

    direc = AltAz(
        location=obs, obstime=Time(time), az=az_deg * u.deg, alt=el_deg * u.deg
    )

    sky = SkyCoord(direc.transform_to(ICRS()))

    return sky.ra.deg, sky.dec.deg


def radec2azel(
    ra_deg: FloatLike,
    dec_deg: FloatLike,
    lat_deg: FloatLike,
    lon_deg: FloatLike,
    time: datetime,
    force_non_astropy: bool = False,
    *,
    delta_ut1: float = 0.0,
) -> tuple:
    """
    sky coordinates (ra, dec) to viewing angle (az, el)

    Parameters
    ----------
    ra_deg : array-like float
         ecliptic right ascension (degress)
    dec_deg : array-like float
         ecliptic declination (degrees)
    lat_deg : array-like float
              observer latitude [-90, 90]
    lon_deg : array-like float
              observer longitude [-180, 180] (degrees)
    time : datetime.datetime or str
           time of observation
    force_non_astropy : bool
        if True, force use of less accurate Numpy implementation even if Astropy is available

    Returns
    -------
    az_deg : array-like float
             azimuth [degrees clockwize from North]
    el_deg : array-like float
             elevation [degrees above horizon (neglecting aberration)]
    """

    if force_non_astropy or "astropy" not in sys.modules:
        return vradec2azel(
            ra_deg, dec_deg, lat_deg, lon_deg, time, delta_ut1=delta_ut1
        )
    else:
        return radec2azel_astropy(ra_deg, dec_deg, lat_deg, lon_deg, time)


def radec2azel_astropy(
    ra_deg: FloatLike,
    dec_deg: FloatLike,
    lat_deg: FloatLike,
    lon_deg: FloatLike,
    time: datetime,
) -> tuple:
    """
    rade2azel using Astropy
    see radec2azel() for description
    """

    obs = EarthLocation(lat=lat_deg * u.deg, lon=lon_deg * u.deg)

    points = SkyCoord(
        Angle(ra_deg, unit=u.deg), Angle(dec_deg, unit=u.deg), equinox="J2000.0"
    )

    altaz = points.transform_to(AltAz(location=obs, obstime=Time(time)))

    return altaz.az.degree, altaz.alt.degree
