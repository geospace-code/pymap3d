"""
converts right ascension, declination to azimuth, elevation and vice versa.
Normally do this via AstroPy.
These functions are fallbacks for those wihtout AstroPy.

Michael Hirsch implementation of algorithms from D. Vallado
"""

import typing
from datetime import datetime

try:
    from numpy import sin, cos, degrees, radians, arcsin as asin, arctan2 as atan2
except ImportError:
    from math import sin, cos, degrees, radians, asin, atan2

from .sidereal import datetime2sidereal

try:
    from numpy.typing import ArrayLike
except ImportError:
    ArrayLike = typing.Any

__all__ = ["azel2radec", "radec2azel"]


def azel2radec(
    az_deg: ArrayLike, el_deg: ArrayLike, lat_deg: ArrayLike, lon_deg: ArrayLike, time: datetime, *, use_astropy: bool = True
) -> typing.Tuple[ArrayLike, ArrayLike]:
    """
    converts azimuth, elevation to right ascension, declination

    Parameters
    ----------

    az_deg : ArrayLike
        azimuth (clockwise) to point [degrees]
    el_deg : ArrayLike
        elevation above horizon to point [degrees]
    lat_deg : ArrayLike
        observer WGS84 latitude [degrees]
    lon_deg : ArrayLike
        observer WGS84 longitude [degrees]
    time : datetime.datetime
        time of observation
    use_astropy : bool, optional
        use AstroPy

    Results
    -------

    ra_deg : ArrayLike
        right ascension to target [degrees]
    dec_deg : ArrayLike
        declination of target [degrees]

    from D.Vallado Fundamentals of Astrodynamics and Applications
    p.258-259
    """

    if abs(lat_deg) > 90:
        raise ValueError("-90 <= lat <= 90")

    az = radians(az_deg)
    el = radians(el_deg)
    lat = radians(lat_deg)
    lon = radians(lon_deg)
    # %% Vallado "algorithm 28" p 268
    dec = asin(sin(el) * sin(lat) + cos(el) * cos(lat) * cos(az))

    lha = atan2(-(sin(az) * cos(el)) / cos(dec), (sin(el) - sin(lat) * sin(dec)) / (cos(dec) * cos(lat)))

    lst = datetime2sidereal(time, lon, use_astropy=use_astropy)  # lon, ra in RADIANS

    """ by definition right ascension [0, 360) degrees """
    return degrees(lst - lha) % 360, degrees(dec)


def radec2azel(
    ra_deg: ArrayLike, dec_deg: ArrayLike, lat_deg: ArrayLike, lon_deg: ArrayLike, time: datetime, *, use_astropy: bool = True
) -> typing.Tuple[ArrayLike, ArrayLike]:
    """
    converts right ascension, declination to azimuth, elevation

    Parameters
    ----------

    ra_deg : ArrayLike
        right ascension to target [degrees]
    dec_deg : ArrayLike
        declination to target [degrees]
    lat_deg : ArrayLike
        observer WGS84 latitude [degrees]
    lon_deg : ArrayLike
        observer WGS84 longitude [degrees]
    time : datetime.datetime
        time of observation
    use_astropy : bool, optional
        use Astropy if available

    Results
    -------

    az_deg : ArrayLike
        azimuth clockwise from north to point [degrees]
    el_deg : ArrayLike
        elevation above horizon to point [degrees]


    from D. Vallado "Fundamentals of Astrodynamics and Applications "
       4th Edition Ch. 4.4 pg. 266-268
    """
    if abs(lat_deg) > 90:
        raise ValueError("-90 <= lat <= 90")

    ra = radians(ra_deg)
    dec = radians(dec_deg)
    lat = radians(lat_deg)
    lon = radians(lon_deg)

    lst = datetime2sidereal(time, lon, use_astropy=use_astropy)  # RADIANS
    # %% Eq. 4-11 p. 267 LOCAL HOUR ANGLE
    lha = lst - ra
    # %% #Eq. 4-12 p. 267
    el = asin(sin(lat) * sin(dec) + cos(lat) * cos(dec) * cos(lha))
    # %% combine Eq. 4-13 and 4-14 p. 268
    az = atan2(-sin(lha) * cos(dec) / cos(el), (sin(dec) - sin(el) * sin(lat)) / (cos(el) * cos(lat)))

    return degrees(az) % 360.0, degrees(el)
