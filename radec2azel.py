#!/usr/bin/env python
"""
Example Kitt Peak

./radec2azel.py 257.96295344 15.43785495 31.9583 -111.5967 2014-12-25T22:00:00MST
"""
from pymap3d import radec2azel

if __name__ == '__main__': 
    from argparse import ArgumentParser
    p = ArgumentParser(description="convert RightAscension,Declination to Azimuth,Elevation")
    p.add_argument('ra',help='right ascension [degrees]',type=float)
    p.add_argument('dec',help='declination [degrees]',type=float)
    p.add_argument('lat',help='WGS84 latitude of observer [degrees]',type=float)
    p.add_argument('lon',help='WGS84 latitude of observer [degrees]',type=float)
    p.add_argument('time',help='UTC time of observation YYYY-mm-ddTHH:MM:SSZ')
    p = p.parse_args()

    az_deg, el_deg = radec2azel(p.ra, p.dec, p.lat, p.lon, p.time)
    print('azimuth: [deg]',az_deg)
    print('elevation [deg]:',el_deg)
