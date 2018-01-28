#!/usr/bin/env python

# Copyright (c) 2014-2018 Michael Hirsch, Ph.D.
#
# Redistribution and use in source and binary forms, with or without modification, are permitted provided that the following conditions are met:
# 1. Redistributions of source code must retain the above copyright notice, this list of conditions and the following disclaimer.
# 2. Redistributions in binary form must reproduce the above copyright notice, this list of conditions and the following disclaimer in the documentation and/or other materials provided with the distribution.
# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.


"""
Example Kitt Peak

./demo_azel2radec.py 264.91833333333335 37.91138888888889 31.9583 -111.5967 2014-12-25T22:00:00MST
"""
from pymap3d import azel2radec

if __name__ == "__main__":  # pragma: no cover
    from argparse import ArgumentParser

    p = ArgumentParser(description="convert azimuth and elevation to "
                       "right ascension and declination")
    p.add_argument('azimuth', help='azimuth [deg]', type=float)
    p.add_argument('elevation', help='elevation [deg]', type=float)
    p.add_argument('lat', help='WGS84 obs. lat [deg]', type=float)
    p.add_argument('lon', help='WGS84 obs. lon [deg]',type=float)
    p.add_argument('time', help='obs. time YYYY-mm-ddTHH:MM:SSZ')
    p = p.parse_args()

    ra, dec = azel2radec(p.azimuth, p.elevation, p.lat, p.lon, p.time)

    print('ra [deg] {} , dec [deg] = {}'.format(ra, dec))
