#!/usr/bin/env python

# Copyright (c) 2014-2018 Michael Hirsch, Ph.D.
#
# Redistribution and use in source and binary forms, with or without modification, are permitted provided that the following conditions are met:
# 1. Redistributions of source code must retain the above copyright notice, this list of conditions and the following disclaimer.
# 2. Redistributions in binary form must reproduce the above copyright notice, this list of conditions and the following disclaimer in the documentation and/or other materials provided with the distribution.
# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

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
