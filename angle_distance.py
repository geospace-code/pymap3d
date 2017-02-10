#!/usr/bin/env python
from numpy.testing import assert_allclose
from pymap3d.haversine import angledist, angledist_meeus


if __name__ == '__main__':  # pragma: no cover
    from argparse import ArgumentParser
    p = ArgumentParser(description="computes angular distance between two points in sky")
    p.add_argument('r0', help='right ascension: first point [deg]', type=float)
    p.add_argument('d0', help='declination: first point [deg]', type=float)
    p.add_argument('r1', help='right ascension: 2nd point [deg]', type=float)
    p.add_argument('d1', help='declination: 2nd point [degrees]', type=float)
    a = p.parse_args()

    dist_deg = angledist_meeus(a.r0, a.d0, a.r1, a.d1)
    dist_deg_astropy = angledist(a.r0, a.d0, a.r1, a.d1)

    print('{:.6f} deg sep'.format(dist_deg))

    assert_allclose(dist_deg, dist_deg_astropy)