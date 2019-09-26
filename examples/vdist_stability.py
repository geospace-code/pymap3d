#!/usr/bin/env python
from pymap3d.vincenty import vdist
import sys
from math import isclose, nan
import typing

eng = None  # don't start engine over and over when script is interactive
try:
    import matlab.engine

    if eng is None:
        eng = matlab.engine.start_matlab("-nojvm")
except Exception as exc:
    print(exc, file=sys.stderr)
    eng = None


def matlab_func(lat1: float, lon1: float, lat2: float, lon2: float) -> typing.Tuple[float, float]:
    """ Using Matlab Engine to do same thing as Pymap3d """
    ell = eng.wgs84Ellipsoid()
    return eng.distance(lat1, lon1, lat2, lon2, ell, nargout=2)


dlast, alast = nan, nan
lon1, lon2 = 0.0, 1.0
for i in range(20):
    lat1 = lat2 = 10.0 ** (-i)
    try:
        dist_m, az_deg = vdist(lat1, lon1, lat2, lon2)
    except Exception as exc:
        print(exc, f"at latitudes {lat1} {lat2}", file=sys.stderr)
        continue
    assert dist_m != dlast
    assert az_deg != alast
    mat_match = True
    dist_matlab, az_matlab = matlab_func(lat1, lon1, lat2, lon2)
    if not isclose(dist_matlab, dist_m):
        mat_match = False
        print(f"MISMATCH: latitude {lat1} {lat2}: Python: {dist_m}  Matlab: {dist_matlab}", file=sys.stderr)
    if not isclose(az_matlab, az_deg):
        mat_match = False
        print(f"MISMATCH: latitude {lat1} {lat2}: Python: {az_matlab}  Matlab: {az_deg}", file=sys.stderr)
    if mat_match:
        print(f"latitudes {lat1} {lat2}: {dist_m} meters {az_deg} deg azimuth")
