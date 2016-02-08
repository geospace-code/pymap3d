#!/usr/bin/env python3
"""
requires Astropy 1.0 or newer. Try usevallado=True if you have an old
Astropy (Vallado accuracy is worse).
Michael Hirsch
GPLv3+
"""
from __future__ import division,absolute_import
import logging
from dateutil.parser import parse
from numpy import sin, cos, degrees, radians,arcsin, arctan2, atleast_1d
from astropy import units as u
from astropy.time import Time
from astropy.coordinates import SkyCoord, EarthLocation, AltAz, ICRS

#from astrometry_azel.datetime2hourangle import datetime2sidereal


def azel2radec(az_deg, el_deg, lat_deg, lon_deg, dtime):
    obs = EarthLocation(lat=lat_deg*u.deg, lon=lon_deg*u.deg)
    direc = AltAz(location=obs, obstime=Time(dtime),
                  az=az_deg*u.deg, alt=el_deg*u.deg)
    sky = SkyCoord(direc.transform_to(ICRS()))

    return sky.ra.deg, sky.dec.deg

def azel2radecvallado(az_deg,el_deg,lat_deg,lon_deg,dtimen):
    """ from D.Vallado Fundamentals of Astrodynamics and Applications
    p.258-259 """
    az_deg = atleast_1d(az_deg)
    el_deg = atleast_1d(el_deg)
    lat_deg = atleast_1d(lat_deg)
    lon_deg = atleast_1d(lon_deg)

    if az_deg.shape != el_deg.shape:
        raise TypeError('az and el must be same shape ndarray')
    if lat_deg.size != 1 or lon_deg.size !=1:
        raise TypeError('need one observer and one or more  (az,el).')

    az = radians(az_deg); el = radians(el_deg)
    lat = radians(lat_deg); lon = radians(lon_deg)
     #Vallado "algorithm 28" p 268
    dec = arcsin( sin(el) * sin(lat) + cos(el) * cos(lat) * cos(az) )

    lha = arctan2( -(sin(az) * cos(el)) / cos(dec),
                   (sin(el) - sin(lat)*sin(dec)) / (cos(dec) * cos(lat)) )

    lst = datetime2sidereal(dtime,lon) #lon, ra in RADIANS

    """ by definition right ascension \in [0,360) degrees """
    return degrees(lst - lha) % 360, degrees(dec)

if __name__ == "__main__": #pragma: no cover
    from argparse import ArgumentParser

    p = ArgumentParser(description="convert azimuth and elevation to "
                 "right ascension and declination")
    p.add_argument('azimuth',help='azimuth [degrees]', nargs='?', type=float)
    p.add_argument('elevation',help='elevation [degrees]', nargs='?', type=float)
    p.add_argument('lat',help='WGS84 latitude of observer [deg] ',nargs='?', type=float)
    p.add_argument('lon',help='WGS84 longitude of observer [deg.]', nargs='?',type=float)
    p.add_argument('time',help='time of observation YYYY-mm-ddTHH:MM:SSZ', nargs='?')
    p.add_argument('--idl',help='run Kitts Peak example from IDL Astrolib',action='store_true')
    a = p.parse_args()



    if a.idl or a.time is None:
        el=37+54/60+41/3600; az=264+55/60+6/3600
        lat=31.9583; lon=-111.5967
        dtime=parse('2014-12-25T22:00:00MST')

        print('demo mode: Kitt Peak.  az {}  el {}  {}'.format(az,el,dtime))
    else:
        dtime = parse(a.time)
        print('using time {}'.format(dtime))

        az=a.azimuth; el=a.elevation; lat=a.lat; lon=a.lon

    ra,dec = azel2radec(az,el,lat,lon,dtime)
    print('ra [deg] {} , dec [deg] = {}'.format(ra,dec) )

