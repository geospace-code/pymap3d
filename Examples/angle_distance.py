#!/usr/bin/env python
from argparse import ArgumentParser

from pymap3d.haversine import anglesep, anglesep_meeus
from pytest import approx

p = ArgumentParser(description="angular distance between two sky points")
p.add_argument("r0", help="right ascension: first point [deg]", type=float)
p.add_argument("d0", help="declination: first point [deg]", type=float)
p.add_argument("r1", help="right ascension: 2nd point [deg]", type=float)
p.add_argument("d1", help="declination: 2nd point [degrees]", type=float)
a = p.parse_args()

dist_deg = anglesep_meeus(a.r0, a.d0, a.r1, a.d1)
dist_deg_astropy = anglesep(a.r0, a.d0, a.r1, a.d1)

print(f"{dist_deg:.6f} deg sep")

assert dist_deg == approx(dist_deg_astropy)
