#!/usr/bin/env python
# Copyright (c) 2014-2018 Michael Hirsch, Ph.D.

"""
converts right ascension, declination to azimuth, elevation and vice versa.
Normally do this via AstroPy.
These functions are fallbacks for those wihtout AstroPy.

Michael Hirsch implementation of algorithms from D. Vallado
"""
from datetime import datetime
from numpy import sin, cos, degrees, radians, arcsin, arctan2, atleast_1d
from typing import Tuple
from .sidereal import datetime2sidereal

__all__ = ['azel2radec', 'radec2azel']


def azel2radec(az_deg: float, el_deg: float,
               lat_deg: float, lon_deg: float,
               time: datetime) -> Tuple[float, float]:
    """
    converts azimuth, elevation to right ascension, declination

    Parameters
    ----------

    az_deg : float or numpy.ndarray of float
        azimuth (clockwise) to point [degrees]

    el_deg : float or numpy.ndarray of float
        elevation above horizon to point [degrees]

    lat_deg : float
        observer WGS84 latitude [degrees]

    lon_deg : float
        observer WGS84 longitude [degrees]

    time : datetime.datetime
        time of observation

    Results
    -------

    ra_deg : float or numpy.ndarray of float
        right ascension to target [degrees]

    dec_deg : float or numpy.ndarray of float
        declination of target [degrees]

    from D.Vallado Fundamentals of Astrodynamics and Applications
    p.258-259
    """
    az = atleast_1d(az_deg)
    el = atleast_1d(el_deg)
    lat = atleast_1d(lat_deg)
    lon = atleast_1d(lon_deg)

    if az.shape != el.shape:
        raise ValueError('az and el must be same shape ndarray')
    if not(lat.size == 1 and lon.size == 1):
        raise ValueError('need one observer and one or more  (az,el).')
    if ((lat < -90) | (lat > 90)).any():
        raise ValueError('-90 <= lat <= 90')

    az = radians(az)
    el = radians(el)
    lat = radians(lat)
    lon = radians(lon)
# %% Vallado "algorithm 28" p 268
    dec = arcsin(sin(el) * sin(lat) + cos(el) * cos(lat) * cos(az))

    lha = arctan2(-(sin(az) * cos(el)) / cos(dec),
                  (sin(el) - sin(lat) * sin(dec)) / (cos(dec) * cos(lat)))

    lst = datetime2sidereal(time, lon)  # lon, ra in RADIANS

    """ by definition right ascension [0, 360) degrees """
    return degrees(lst - lha) % 360, degrees(dec)


def radec2azel(ra_deg: float, dec_deg: float,
               lat_deg: float, lon_deg: float,
               time: datetime) -> Tuple[float, float]:
    """
    converts right ascension, declination to azimuth, elevation

    Parameters
    ----------

    ra_deg : float or numpy.ndarray of float
        right ascension to target [degrees]

    dec_deg : float or numpy.ndarray of float
        declination to target [degrees]

    lat_deg : float
        observer WGS84 latitude [degrees]

    lon_deg : float
        observer WGS84 longitude [degrees]

    time : datetime.datetime
        time of observation

    Results
    -------

    az_deg : float or numpy.ndarray of float
        azimuth clockwise from north to point [degrees]

    el_deg : float or numpy.ndarray of float
        elevation above horizon to point [degrees]


    from D. Vallado "Fundamentals of Astrodynamics and Applications "
       4th Edition Ch. 4.4 pg. 266-268
    """
    ra = atleast_1d(ra_deg)
    dec = atleast_1d(dec_deg)
    lat = atleast_1d(lat_deg)
    lon = atleast_1d(lon_deg)

    if ra.shape != dec.shape:
        raise ValueError('az and el must be same shape ndarray')
    if not(lat.size == 1 and lon.size == 1):
        raise ValueError('need one observer and one or more  (az,el).')
    if ((lat < -90) | (lat > 90)).any():
        raise ValueError('-90 <= lat <= 90')

    ra = radians(ra)
    dec = radians(dec)
    lat = radians(lat)
    lon = radians(lon)

    lst = datetime2sidereal(time, lon)  # RADIANS
# %% Eq. 4-11 p. 267 LOCAL HOUR ANGLE
    lha = lst - ra
# %% #Eq. 4-12 p. 267
    el = arcsin(sin(lat) * sin(dec) + cos(lat) * cos(dec) * cos(lha))
# %% combine Eq. 4-13 and 4-14 p. 268
    az = arctan2(-sin(lha) * cos(dec) / cos(el),
                 (sin(dec) - sin(el) * sin(lat)) / (cos(el) * cos(lat)))

    return degrees(az) % 360.0, degrees(el)
