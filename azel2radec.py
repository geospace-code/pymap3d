#!/usr/bin/env python3
from __future__ import division
from numpy import sin, cos, degrees, radians,arcsin, arctan2,atleast_1d, nan
import sys
sys.path.append('../astrometry') #https://github.com/scienceopen/astrometry/
from datetime2hourangle import datetime2sidereal
#from pdb import set_trace

def azel2radec(az_deg,el_deg,lat_deg,lon_deg,dtime):

    """ from D.Vallado Fundamentals of Astrodynamics and Applications p.258-259 """
    az = atleast_1d(radians(az_deg))
    el = atleast_1d(radians(el_deg))
    lat = atleast_1d(radians(lat_deg))
    lon = atleast_1d(radians(lon_deg))

    if az.shape != el.shape: exit('az and el must be same shape')
    if lat.shape != lat.shape: exit('lat and lon must be same shape')

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
    p.add_argument('--selftest',help='integration test',action='store_true')
    a = p.parse_args()

    if a.selftest:
        from numpy.testing import assert_almost_equal
        ra,dec = azel2radec(180.1, 80, 65, -148,parse('2014-04-06T08:00:00Z'))
        assert_almost_equal(ra,166.5032081149338)
        assert_almost_equal(dec,55.000011165405752)
        exit(0)

    dtime = parse(a.time)
    print(dtime)

    ra,dec = azel2radec(a.azimuth,a.elevation,a.lat,a.lon,dtime)
    print('ra / dec =',(ra,dec) )
