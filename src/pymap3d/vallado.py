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

__all__ = ["azel2radec", "radec2azel"]

if typing.TYPE_CHECKING:
    from numpy import ndarray


def azel2radec(
    az_deg: "ndarray", el_deg: "ndarray", lat_deg: "ndarray", lon_deg: "ndarray", time: datetime, *, use_astropy: bool = True
) -> typing.Tuple["ndarray", "ndarray"]:
    """
    converts azimuth, elevation to right ascension, declination

    Parameters
    ----------

    az_deg : "ndarray"
        azimuth (clockwise) to point [degrees]
    el_deg : "ndarray"
        elevation above horizon to point [degrees]
    lat_deg : "ndarray"
        observer WGS84 latitude [degrees]
    lon_deg : "ndarray"
        observer WGS84 longitude [degrees]
    time : datetime.datetime
        time of observation
    use_astropy : bool, optional
        use AstroPy

    Results
    -------

    ra_deg : "ndarray"
        right ascension to target [degrees]
    dec_deg : "ndarray"
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
    ra_deg: "ndarray", dec_deg: "ndarray", lat_deg: "ndarray", lon_deg: "ndarray", time: datetime, *, use_astropy: bool = True
) -> typing.Tuple["ndarray", "ndarray"]:
    """
    converts right ascension, declination to azimuth, elevation

    Parameters
    ----------

    ra_deg : "ndarray"
        right ascension to target [degrees]
    dec_deg : "ndarray"
        declination to target [degrees]
    lat_deg : "ndarray"
        observer WGS84 latitude [degrees]
    lon_deg : "ndarray"
        observer WGS84 longitude [degrees]
    time : datetime.datetime
        time of observation
    use_astropy : bool, optional
        use Astropy if available

    Results
    -------

    az_deg : "ndarray"
        azimuth clockwise from north to point [degrees]
    el_deg : "ndarray"
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
