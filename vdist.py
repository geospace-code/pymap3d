#!/usr/bin/env python
from pymap3d.vincenty import vdist
from argparse import ArgumentParser


def main():
    p = ArgumentParser(description='vdist distance between WGS-84 coordinates')
    p.add_argument('lat1', help='latitude1 WGS-84 [degrees]', type=float)
    p.add_argument('lon1', help='longitude1 WGS-84 [degrees]', type=float)
    p.add_argument('lat2', help='latitude2 WGS-84 [degrees]', type=float)
    p.add_argument('lon2', help='longitude2 WGS-84 [degrees]', type=float)
    P = p.parse_args()

    dist_m = vdist(P.lat1, P.lon1, P.lat2, P.lon2)
    print('distance between WGS-84 points: {} m'.format(dist_m))


if __name__ == '__main__':  # pragma: no cover
    main()
