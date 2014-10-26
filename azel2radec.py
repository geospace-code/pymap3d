#!/usr/bin/env python3

from numpy import sin, cos, degrees, radians,arcsin, arctan2,atleast_1d
import sys
sys.path.append('../astrometry') #https://github.com/scienceopen/astrometry/
from datetime2hourangle import datetime2sidereal
#from pdb import set_trace

def azel2radec(az_deg,el_deg,lat_deg,lon_deg,dtime):

#from D.Vallado Fundamentals of Astrodynamics and Applications p.258-259
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
    ra  = lst - lha

    return degrees(ra), degrees(dec)

if __name__ == "__main__":
    from dateutil.parser import parse
    from argparse import ArgumentParser

    p = ArgumentParser(description='convert azimuth and elevation to right ascension and declination')
    p.add_argument('azimuth',help='azimuth [degrees]',type=float)
    p.add_argument('elevation',help='elevation [degrees]',type=float)
    p.add_argument('latitude',help='WGS84 latitude of observer [deg] ',type=float)
    p.add_argument('longitude',help='WGS84 longitude of observer [deg.]',type=float)
    p.add_argument('time',help='time of observation YYYY-mm-ddTHH:MM:SSZ',type=str)
    args = p.parse_args()

    az = args.azimuth
    el = args.elevation
    lat = args.latitude
    lon = args.longitude
    time = args.time

    dtime = parse(time)

    print(dtime)

    ra,dec = azel2radec(az,el,lat,lon,dtime)


    print('ra / dec = ' + str((ra,dec)) )