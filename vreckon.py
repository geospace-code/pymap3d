#!/usr/bin/env python
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
    print('new lat / lon {} / {} '.format(lat2,lon2))
    print('az back to start:', a21)