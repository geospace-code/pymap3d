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
from .datetime2hourangle import datetime2sidereal


def azel2radec(az_deg: float, el_deg: float,
               lat_deg: float, lon_deg: float,
               time: datetime) -> Tuple[float, float]:
    """
    convert azimuth, elevation to right ascension, declination

    Inputs

    az_deg
        Numpy ndarray of azimuth to point [degrees]

    el_deg
        Numpy ndarray of elevation to point [degrees]

    lat_deg
        scalar observer WGS84 latitude [degrees]

    lon_deg
        scalar observer WGS84 longitude [degrees]

    time
        time of observation

    Outputs

    ra_deg
        Numpy ndarray of right ascension values [degrees]

    dec_deg
        Numpy ndarray of declination values [degrees]

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
    if ((lon < -180) | (lon > 360)).any():
        raise ValueError('-180 <= lat <= 360')

    az = radians(az)
    el = radians(el)
    lat = radians(lat)
    lon = radians(lon)
# %% Vallado "algorithm 28" p 268
    dec = arcsin(sin(el) * sin(lat) + cos(el) * cos(lat) * cos(az))

    lha = arctan2(-(sin(az) * cos(el)) / cos(dec),
                  (sin(el) - sin(lat) * sin(dec)) / (cos(dec) * cos(lat)))

    lst = datetime2sidereal(time, lon)  # lon, ra in RADIANS

    """ by definition right ascension \in [0,360) degrees """
    return degrees(lst - lha) % 360, degrees(dec)


def radec2azel(ra_deg: float, dec_deg: float,
               lat_deg: float, lon_deg: float,
               time: datetime) -> Tuple[float, float]:
    """
    convert right ascension, declination to azimuth, elevation

    Inputs

    ra_deg
        Numpy ndarray of right ascension values [degrees]

    dec_deg
        Numpy ndarray of declination values [degrees]

    lat_deg
        scalar observer WGS84 latitude [degrees]

    lon_deg
        scalar observer WGS84 longitude [degrees]

    time
        time of observation

    Outputs

    az_deg
        Numpy ndarray of azimuth to point [degrees]

    el_deg
        Numpy ndarray of elevation to point [degrees]


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
    if ((lon < -180) | (lon > 360)).any():
        raise ValueError('-180 <= lat <= 360')

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
