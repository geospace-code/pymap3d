#!/usr/bin/env python

# Copyright (c) 2014-2018 Michael Hirsch, Ph.D.
#
# Redistribution and use in source and binary forms, with or without modification, are permitted provided that the following conditions are met:
# 1. Redistributions of source code must retain the above copyright notice, this list of conditions and the following disclaimer.
# 2. Redistributions in binary form must reproduce the above copyright notice, this list of conditions and the following disclaimer in the documentation and/or other materials provided with the distribution.
# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

from pymap3d.vincenty import vdist

if __name__ == '__main__': #pragma: no cover
    from argparse import ArgumentParser
    p = ArgumentParser(description='vdist distance between WGS-84 coordinates')
    p.add_argument('lat1',help='latitude1 WGS-84 [degrees]',type=float)
    p.add_argument('lon1',help='longitude1 WGS-84 [degrees]',type=float)
    p.add_argument('lat2',help='latitude2 WGS-84 [degrees]',type=float)
    p.add_argument('lon2',help='longitude2 WGS-84 [degrees]',type=float)
    p = p.parse_args()

    dist_m = vdist(p.lat1, p.lon1, p.lat2, p.lon2)
    print('distance between WGS-84 points: {} m'.format(dist_m))
