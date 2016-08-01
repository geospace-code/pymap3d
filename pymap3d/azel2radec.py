#!/usr/bin/env python
"""
requires Astropy 1.0 or newer. Try usevallado=True if you have an old
Astropy (Vallado accuracy is worse).
"""
from __future__ import division
from six import string_types
from numpy import sin, cos, degrees, radians, arcsin, arctan2, atleast_1d
from astropy import units as u
from astropy.time import Time
from astropy.coordinates import SkyCoord, EarthLocation, AltAz, ICRS
#
from .datetime2hourangle import datetime2sidereal
from .common import str2dt

def azel2radec(az_deg, el_deg, lat_deg, lon_deg, t):

    t = str2dt(t)

    obs = EarthLocation(lat=lat_deg * u.deg, lon=lon_deg * u.deg)
    direc = AltAz(location=obs, obstime=Time(t),
                  az=az_deg * u.deg, alt=el_deg * u.deg)
    sky = SkyCoord(direc.transform_to(ICRS()))

    return sky.ra.deg, sky.dec.deg


def azel2radecvallado(az_deg, el_deg, lat_deg, lon_deg, t):
    """ from D.Vallado Fundamentals of Astrodynamics and Applications
    p.258-259 """
    az_deg = atleast_1d(az_deg)
    el_deg = atleast_1d(el_deg)
    lat_deg = atleast_1d(lat_deg)
    lon_deg = atleast_1d(lon_deg)

    assert az_deg.shape == el_deg.shape, 'az and el must be same shape ndarray'
    assert lat_deg.size == 1 and lon_deg.size == 1, 'need one observer and one or more  (az,el).'

    az = radians(az_deg)
    el = radians(el_deg)
    lat = radians(lat_deg)
    lon = radians(lon_deg)
    # Vallado "algorithm 28" p 268
    dec = arcsin(sin(el) * sin(lat) + cos(el) * cos(lat) * cos(az))

    lha = arctan2(-(sin(az) * cos(el)) / cos(dec),
                  (sin(el) - sin(lat) * sin(dec)) / (cos(dec) * cos(lat)))

    lst = datetime2sidereal(t, lon)  # lon, ra in RADIANS

    """ by definition right ascension \in [0,360) degrees """
    return degrees(lst - lha) % 360, degrees(dec)
