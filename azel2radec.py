#!/usr/bin/env python3
"""
requires Astropy 1.0 or newer. Try usevallado=True if you have an old Astropy (Vallado accuracy is worse).
Michael Hirsch
GPLv3+
"""
from __future__ import division
from numpy import sin, cos, degrees, radians,arcsin, arctan2, atleast_1d, nan
import sys


try:
    from astropy import units as u
    from astropy.time import Time
    from astropy.coordinates import SkyCoord, EarthLocation, AltAz, ICRS
    usevallado=False
except ImportError as e:
    print(str(e) + ' trouble importing AstroPy>1.0, falling back to Vallado')
    sys.path.append('../astrometry') # git clone https://github.com/scienceopen/astrometry/
    from datetime2hourangle import datetime2sidereal
    usevallado=True


def azel2radec(az_deg, el_deg, lat_deg, lon_deg, dtime):
    """ from D.Vallado Fundamentals of Astrodynamics and Applications p.258-259 """
    az_deg = atleast_1d(az_deg)
    el_deg = atleast_1d(el_deg)
    lat_deg = atleast_1d(lat_deg)
    lon_deg = atleast_1d(lon_deg)

    if az_deg.shape != el_deg.shape:
        exit('*** azel2radec: az and el must be same shape ndarray')
    if lat_deg.size != 1 or lon_deg.size !=1:
        exit('*** azel2radec is designed for one observer and one or more points (az,el).')

    if usevallado:
        ra_deg, dec_deg = azel2radecvallado(az_deg,el_deg,lat_deg,lon_deg,dtime)
    else: #use astropy v1.0 +
        obs = EarthLocation(lat=lat_deg*u.deg, lon=lon_deg*u.deg)
        direc = AltAz(location=obs, obstime=Time(dtime),
                      az=az_deg*u.deg, alt=el_deg*u.deg)
        sky = SkyCoord(direc.transform_to(ICRS()))

    return sky.ra.deg, sky.dec.deg

def azel2radecvallado(az_deg,el_deg,lat_deg,lon_deg,dtimen):
    az = radians(az_deg); el = radians(el_deg)
    lat = radians(lat_deg); lon = radians(lon_deg)
     #Vallado "algorithm 28" p 268
    dec = arcsin( sin(el) * sin(lat) + cos(el) * cos(lat) * cos(az) )

    lha = arctan2( -(sin(az) * cos(el)) / cos(dec),
                   (sin(el) - sin(lat)*sin(dec)) / (cos(dec) * cos(lat)) )

    lst = datetime2sidereal(dtime,lon) #lon, ra in RADIANS

    """ by definition right ascension \in [0,360) degrees """
    return degrees(lst - lha) % 360, degrees(dec)

if __name__ == "__main__":
    from dateutil.parser import parse
    from argparse import ArgumentParser

    p = ArgumentParser(description='convert azimuth and elevation to right ascension and declination')
    p.add_argument('azimuth',help='azimuth [degrees]',nargs='?',type=float,default=nan)
    p.add_argument('elevation',help='elevation [degrees]',nargs='?',type=float,default=nan)
    p.add_argument('lat',help='WGS84 latitude of observer [deg] ',nargs='?',type=float,default=nan)
    p.add_argument('lon',help='WGS84 longitude of observer [deg.]',nargs='?',type=float,default=nan)
    p.add_argument('time',help='time of observation YYYY-mm-ddTHH:MM:SSZ',nargs='?',type=str,default='')
    a = p.parse_args()

    dtime = parse(a.time)
    print(dtime)

    ra,dec = azel2radec(a.azimuth,a.elevation,a.lat,a.lon,dtime)
    print('ra / dec =',(ra,dec) )
