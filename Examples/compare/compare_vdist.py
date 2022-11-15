#!/usr/bin/env python3
from __future__ import annotations

import sys
from math import isclose, nan
import numpy as np
from pathlib import Path

from matlab_mapping import matlab_mapping

import matlab.engine
from pymap3d.vincenty import vdist

cwd = Path(__file__).parent

eng = matlab.engine.start_matlab("-nojvm")
eng.addpath(eng.genpath(str(cwd)), nargout=0)

has_map = matlab_mapping(eng)


def distance(lat1, lon1, lat2, lon2) -> tuple[float, float]:
    """Using Matlab Engine to do same thing as Pymap3d"""

    if has_map:
        return eng.distance(lat1, lon1, lat2, lon2, eng.wgs84Ellipsoid(), nargout=2)
    else:
        return eng.matmap3d.vdist(lat1, lon1, lat2, lon2, nargout=2)


dlast, alast = nan, nan
lon1, lon2 = 0.0, 1.0
for i in range(20):
    lat1 = lat2 = 10.0 ** (-i)

    dist_m, az_deg = vdist(lat1, lon1, lat2, lon2)

    assert dist_m != dlast
    assert az_deg != alast
    mat_match = True
    dist_matlab, az_matlab = distance(lat1, lon1, lat2, lon2)

    assert np.isreal(dist_m), f"Python vdist distance is not real for input latitude {lat1}"
    assert np.isreal(dist_matlab), f"Matlab distance is not real for input latitude {lat1}"

    if not isclose(dist_matlab, dist_m):
        mat_match = False
        print(
            f"MISMATCH: latitude {lat1} {lat2}: Python: {dist_m}  Matlab: {dist_matlab}",
            file=sys.stderr,
        )
    if not isclose(az_matlab, az_deg):
        mat_match = False
        print(
            f"MISMATCH: latitude {lat1} {lat2}: Python: {az_matlab}  Matlab: {az_deg}",
            file=sys.stderr,
        )
    if mat_match:
        print(f"latitudes {lat1} {lat2}: {dist_m} meters {az_deg} deg azimuth")
