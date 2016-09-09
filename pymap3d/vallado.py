#!/usr/bin/env python
"""
converts right ascension, declination to azimuth, elevation and vice versa



inputs:
ra_deg: numpy ndarray of right ascension values [degrees]
dec_deg: numpy ndarray of declination values [degrees]
lat_deg: scalar observer WGS84 latitude [degrees]
lon_deg: scalar observer WGS84 longitude [degrees]
dtime: UTC time of observation YYYY-mm-ddTHH:MM:SS

outputs:
az_deg: numpy ndarray of azimuth to point [degrees]
el_deg: numpy ndarray of elevation to point [degrees]

Michael Hirsch implementation of algorithms from D. Vallado
"""
from __future__ import division
from numpy import sin, cos, degrees, radians, arcsin, arctan2, atleast_1d

#
from .datetime2hourangle import datetime2sidereal


def azel2radec(az_deg, el_deg, lat_deg, lon_deg, t):
    """
    from D.Vallado Fundamentals of Astrodynamics and Applications
    p.258-259
    """
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


def radec2azel(dtime,ra_deg,dec_deg,lat_deg,lon_deg):
    """
    from D. Vallado "Fundamentals of Astrodynamics and Applications "
       4th Edition Ch. 4.4 pg. 266-268
    """
    ra =  radians(ra_deg)
    dec = radians(dec_deg)
    lat = radians(lat_deg)
    lon = radians(lon_deg)

    lst = datetime2sidereal(dtime,lon) #RADIANS
    lha = lst - ra #Eq. 4-11 p. 267 LOCAL HOUR ANGLE

    el = arcsin(sin(lat) * sin(dec) + cos(lat) * cos(dec) * cos(lha) ) #Eq. 4-12 p. 267
    #combine Eq. 4-13 and 4-14 p. 268
    az = arctan2( -sin(lha) * cos(dec) / cos(el),
                   ( sin(dec) - sin(el) * sin(lat) ) / (cos(el)*cos(lat))  )


    return degrees(az) % 360.0, degrees(el)