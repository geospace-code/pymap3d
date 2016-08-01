#!/usr/bin/env python
"""
A thin shell vdist() around GeoPy Vincenty distance algorithm for those who might be looking
for the vdist() companion to the Michael Kleder vreckon()
"""
from geopy.distance import vincenty

def vdist(lat1,lon1,lat2,lon2):
    return vincenty((lat1,lon1),(lat2,lon2)).meters
