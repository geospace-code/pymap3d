#!/usr/bin/env python3

from numpy import sin, cos, degrees, radians,arcsin, arctan2,atleast_1d
import sys
sys.path.append('../astrometry')
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
    from numpy import linspace
    from dateutil.parser import parse
    #test function
    az = linspace(0,360,10)
    el = linspace(0,90,10)
    lat = 40
    lon = -140
    dtime = parse('2014-08-08T12:13:14Z')

    ra,dec = azel2radec(az,el,lat,lon,dtime)

    print('testing!')
    print('ra / dec = ' + str((ra,dec)) )