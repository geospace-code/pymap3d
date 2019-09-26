#!/usr/bin/env python
from pymap3d.vincenty import vdist


for i in range(1, 16):
    lat = 10**(-i)
    dist_m, az_deg = vdist(lat, 0, -lat, 0)
    print(f"latitude {lat}: {dist_m} meters {az_deg} deg azimuth")
