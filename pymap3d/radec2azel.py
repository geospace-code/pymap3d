#!/usr/bin/env python
"""
converts right ascension, declination to azimuth, elevation

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

New feature with Astropy v1.0-- more accurate calculation, but it's 21 times slower than Vallado.
We're only talking milliseconds here, so I use the more accurate Astropy via usevallado=False
"""
from __future__ import division
from numpy import radians, sin, cos, degrees, arcsin, arctan2, atleast_1d
from astropy import units as u
from astropy.time import Time
from astropy.coordinates import Angle,SkyCoord, EarthLocation, AltAz
#
from .datetime2hourangle import datetime2sidereal
from .common import str2dt

def radec2azel(ra_deg,dec_deg,lat_deg,lon_deg, t):
    """
     for azel2radec please see:
     https://github.com/scienceopen/python-mapping/

    from D. Vallado "Fundamentals of Astrodynamics and Applications "
          4th Edition Ch. 4.4 pg. 266-268
    """
#%% input trapping
    t = str2dt(t)
    lat_deg = atleast_1d(lat_deg); lon_deg = atleast_1d(lon_deg)
    ra_deg = atleast_1d(ra_deg); dec_deg = atleast_1d(dec_deg)

    assert lat_deg.size == 1 & lon_deg.size ==1, 'radec2azel is designed for one observer and one or more points (ra,dec).'
    assert ra_deg.shape == dec_deg.shape, 'ra and dec must be the same shape ndarray'

    obs = EarthLocation(lat=lat_deg*u.deg, lon=lon_deg*u.deg)
    points = SkyCoord(Angle(ra_deg, unit=u.deg), Angle(dec_deg, unit=u.deg), equinox='J2000.0')
    altaz = points.transform_to(AltAz(location=obs, obstime=Time(t)))

    return altaz.az.degree, altaz.alt.degree

def radec2azelvallado(dtime,ra_deg,dec_deg,lat_deg,lon_deg):
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
