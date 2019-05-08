"""
Azimuth / elevation <==> Right ascension, declination
"""
from typing import Tuple
from datetime import datetime
from .vallado import azel2radec as vazel2radec, radec2azel as vradec2azel
from .timeconv import str2dt  # astropy can't handle xarray times (yet)
try:
    from astropy.time import Time
    from astropy import units as u
    from astropy.coordinates import Angle, SkyCoord, EarthLocation, AltAz, ICRS
except ImportError:
    Time = None


def azel2radec(az_deg: float, el_deg: float,
               lat_deg: float, lon_deg: float,
               time: datetime, usevallado: bool = False) -> Tuple[float, float]:
    """
    viewing angle (az, el) to sky coordinates (ra, dec)

    Parameters
    ----------
    az_deg : float or numpy.ndarray of float
         azimuth [degrees clockwize from North]
    el_deg : float or numpy.ndarray of float
             elevation [degrees above horizon (neglecting aberration)]
    lat_deg : float
              observer latitude [-90, 90]
    lon_deg : float
              observer longitude [-180, 180] (degrees)
    time : datetime.datetime or str
           time of observation
    usevallado : bool, optional
                 default use astropy. If true, use Vallado algorithm

    Returns
    -------
    ra_deg : float or numpy.ndarray of float
         ecliptic right ascension (degress)
    dec_deg : float or numpy.ndarray of float
         ecliptic declination (degrees)
    """

    if usevallado or Time is None:  # non-AstroPy method, less accurate
        return vazel2radec(az_deg, el_deg, lat_deg, lon_deg, time)

    obs = EarthLocation(lat=lat_deg * u.deg, lon=lon_deg * u.deg)

    direc = AltAz(location=obs, obstime=Time(str2dt(time)),
                  az=az_deg * u.deg, alt=el_deg * u.deg)

    sky = SkyCoord(direc.transform_to(ICRS()))

    return sky.ra.deg, sky.dec.deg


def radec2azel(ra_deg: float, dec_deg: float,
               lat_deg: float, lon_deg: float,
               time: datetime, usevallado: bool = False) -> Tuple[float, float]:
    """
    sky coordinates (ra, dec) to viewing angle (az, el)

    Parameters
    ----------
    ra_deg : float or numpy.ndarray of float
         ecliptic right ascension (degress)
    dec_deg : float or numpy.ndarray of float
         ecliptic declination (degrees)
    lat_deg : float
              observer latitude [-90, 90]
    lon_deg : float
              observer longitude [-180, 180] (degrees)
    time : datetime.datetime or str
           time of observation
    usevallado : bool, optional
                 default use astropy. If true, use Vallado algorithm

    Returns
    -------
    az_deg : float or numpy.ndarray of float
             azimuth [degrees clockwize from North]
    el_deg : float or numpy.ndarray of float
             elevation [degrees above horizon (neglecting aberration)]
    """

    if usevallado or Time is None:
        return vradec2azel(ra_deg, dec_deg, lat_deg, lon_deg, time)

    obs = EarthLocation(lat=lat_deg * u.deg,
                        lon=lon_deg * u.deg)

    points = SkyCoord(Angle(ra_deg, unit=u.deg),
                      Angle(dec_deg, unit=u.deg),
                      equinox='J2000.0')

    altaz = points.transform_to(AltAz(location=obs, obstime=Time(str2dt(time))))

    return altaz.az.degree, altaz.alt.degree
