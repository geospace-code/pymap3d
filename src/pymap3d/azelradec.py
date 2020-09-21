"""
Azimuth / elevation <==> Right ascension, declination
"""
import typing
from datetime import datetime
from .vallado import azel2radec as vazel2radec, radec2azel as vradec2azel
from .timeconv import str2dt  # astropy can't handle xarray times (yet)

try:
    from astropy.time import Time
    from astropy import units as u
    from astropy.coordinates import Angle, SkyCoord, EarthLocation, AltAz, ICRS
except ImportError:
    Time = None

try:
    from numpy.typing import ArrayLike
except ImportError:
    ArrayLike = typing.Any

__all__ = ["radec2azel", "azel2radec"]


def azel2radec(
    az_deg: ArrayLike, el_deg: ArrayLike, lat_deg: ArrayLike, lon_deg: ArrayLike, time: datetime, *, use_astropy: bool = True
) -> typing.Tuple[ArrayLike, ArrayLike]:
    """
    viewing angle (az, el) to sky coordinates (ra, dec)

    Parameters
    ----------
    az_deg : ArrayLike
         azimuth [degrees clockwize from North]
    el_deg : ArrayLike
             elevation [degrees above horizon (neglecting aberration)]
    lat_deg : ArrayLike
              observer latitude [-90, 90]
    lon_deg : ArrayLike
              observer longitude [-180, 180] (degrees)
    time : datetime.datetime or str
           time of observation
    use_astropy : bool, optional
                 default use astropy.

    Returns
    -------
    ra_deg : ArrayLike
         ecliptic right ascension (degress)
    dec_deg : ArrayLike
         ecliptic declination (degrees)
    """

    if use_astropy and Time is not None:

        obs = EarthLocation(lat=lat_deg * u.deg, lon=lon_deg * u.deg)

        direc = AltAz(location=obs, obstime=Time(str2dt(time)), az=az_deg * u.deg, alt=el_deg * u.deg)

        sky = SkyCoord(direc.transform_to(ICRS()))

        return sky.ra.deg, sky.dec.deg

    return vazel2radec(az_deg, el_deg, lat_deg, lon_deg, time)


def radec2azel(
    ra_deg: ArrayLike, dec_deg: ArrayLike, lat_deg: ArrayLike, lon_deg: ArrayLike, time: datetime, *, use_astropy: bool = False
) -> typing.Tuple[ArrayLike, ArrayLike]:
    """
    sky coordinates (ra, dec) to viewing angle (az, el)

    Parameters
    ----------
    ra_deg : ArrayLike
         ecliptic right ascension (degress)
    dec_deg : ArrayLike
         ecliptic declination (degrees)
    lat_deg : ArrayLike
              observer latitude [-90, 90]
    lon_deg : ArrayLike
              observer longitude [-180, 180] (degrees)
    time : datetime.datetime or str
           time of observation
    use_astropy : bool, optional
                 default use astropy.

    Returns
    -------
    az_deg : ArrayLike
             azimuth [degrees clockwize from North]
    el_deg : ArrayLike
             elevation [degrees above horizon (neglecting aberration)]
    """

    if use_astropy and Time is not None:
        obs = EarthLocation(lat=lat_deg * u.deg, lon=lon_deg * u.deg)
        points = SkyCoord(Angle(ra_deg, unit=u.deg), Angle(dec_deg, unit=u.deg), equinox="J2000.0")
        altaz = points.transform_to(AltAz(location=obs, obstime=Time(str2dt(time))))

        return altaz.az.degree, altaz.alt.degree

    return vradec2azel(ra_deg, dec_deg, lat_deg, lon_deg, time)
