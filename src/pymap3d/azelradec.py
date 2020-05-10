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

__all__ = ["radec2azel", "azel2radec"]

if typing.TYPE_CHECKING:
    from numpy import ndarray


def azel2radec(
    az_deg: "ndarray", el_deg: "ndarray", lat_deg: "ndarray", lon_deg: "ndarray", time: datetime, *, use_astropy: bool = True
) -> typing.Tuple["ndarray", "ndarray"]:
    """
    viewing angle (az, el) to sky coordinates (ra, dec)

    Parameters
    ----------
    az_deg : "ndarray"
         azimuth [degrees clockwize from North]
    el_deg : "ndarray"
             elevation [degrees above horizon (neglecting aberration)]
    lat_deg : "ndarray"
              observer latitude [-90, 90]
    lon_deg : "ndarray"
              observer longitude [-180, 180] (degrees)
    time : datetime.datetime or str
           time of observation
    use_astropy : bool, optional
                 default use astropy.

    Returns
    -------
    ra_deg : "ndarray"
         ecliptic right ascension (degress)
    dec_deg : "ndarray"
         ecliptic declination (degrees)
    """

    if use_astropy and Time is not None:

        obs = EarthLocation(lat=lat_deg * u.deg, lon=lon_deg * u.deg)

        direc = AltAz(location=obs, obstime=Time(str2dt(time)), az=az_deg * u.deg, alt=el_deg * u.deg)

        sky = SkyCoord(direc.transform_to(ICRS()))

        return sky.ra.deg, sky.dec.deg

    return vazel2radec(az_deg, el_deg, lat_deg, lon_deg, time)


def radec2azel(
    ra_deg: "ndarray", dec_deg: "ndarray", lat_deg: "ndarray", lon_deg: "ndarray", time: datetime, *, use_astropy: bool = False
) -> typing.Tuple["ndarray", "ndarray"]:
    """
    sky coordinates (ra, dec) to viewing angle (az, el)

    Parameters
    ----------
    ra_deg : "ndarray"
         ecliptic right ascension (degress)
    dec_deg : "ndarray"
         ecliptic declination (degrees)
    lat_deg : "ndarray"
              observer latitude [-90, 90]
    lon_deg : "ndarray"
              observer longitude [-180, 180] (degrees)
    time : datetime.datetime or str
           time of observation
    use_astropy : bool, optional
                 default use astropy.

    Returns
    -------
    az_deg : "ndarray"
             azimuth [degrees clockwize from North]
    el_deg : "ndarray"
             elevation [degrees above horizon (neglecting aberration)]
    """

    if use_astropy and Time is not None:
        obs = EarthLocation(lat=lat_deg * u.deg, lon=lon_deg * u.deg)
        points = SkyCoord(Angle(ra_deg, unit=u.deg), Angle(dec_deg, unit=u.deg), equinox="J2000.0")
        altaz = points.transform_to(AltAz(location=obs, obstime=Time(str2dt(time))))

        return altaz.az.degree, altaz.alt.degree

    return vradec2azel(ra_deg, dec_deg, lat_deg, lon_deg, time)
