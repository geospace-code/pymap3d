#!/usr/bin/env python
"""
Example Kitt Peak

./demo_azel2radec.py 264.9183 37.911388 31.9583 -111.597 2014-12-25T22:00:00MST
"""
from argparse import ArgumentParser

from pymap3d import azel2radec

p = ArgumentParser(
    description="convert azimuth and elevation to " "right ascension and declination"
)
p.add_argument("azimuth", help="azimuth [deg]", type=float)
p.add_argument("elevation", help="elevation [deg]", type=float)
p.add_argument("lat", help="WGS84 obs. lat [deg]", type=float)
p.add_argument("lon", help="WGS84 obs. lon [deg]", type=float)
p.add_argument("time", help="obs. time YYYY-mm-ddTHH:MM:SSZ")
P = p.parse_args()

ra, dec = azel2radec(P.azimuth, P.elevation, P.lat, P.lon, P.time)

print("ra [deg] ", ra, " dec [deg] ", dec)
