#!/usr/bin/env python

# Copyright (c) 2014-2018 Michael Hirsch, Ph.D.
#
# Redistribution and use in source and binary forms, with or without modification, are permitted provided that the following conditions are met:
# 1. Redistributions of source code must retain the above copyright notice, this list of conditions and the following disclaimer.
# 2. Redistributions in binary form must reproduce the above copyright notice, this list of conditions and the following disclaimer in the documentation and/or other materials provided with the distribution.
# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

from pymap3d.vincenty import vreckon

if __name__ == '__main__':  # pragma: no cover
    from argparse import ArgumentParser

    p = ArgumentParser(description='Python port of vreckon.m')
    p.add_argument('lat', help='latitude WGS-84 [degrees]', type=float)
    p.add_argument('lon', help='longitude WGS-84 [degrees]', type=float)
    p.add_argument('range', help='range from start point [meters]', type=float)
    p.add_argument('azimuth', help='azimuth to start [deg.]', type=float)
    p = p.parse_args()

    lat2, lon2, a21 = vreckon(p.lat, p.lon, p.range, p.azimuth)
    print('new (lat, lon)  ({}, {}) '.format(lat2,lon2))
    print('az back to start:', a21)
