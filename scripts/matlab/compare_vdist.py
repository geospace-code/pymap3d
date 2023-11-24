#!/usr/bin/env python3
"""
Compare with Matlab Mapping Toolbox distance()
"""

from __future__ import annotations

import logging
from math import isclose, nan
import numpy as np

from matlab_engine import matlab_engine, has_mapping

from pymap3d.vincenty import vdist


def distance(eng, lat1, lon1, lat2, lon2) -> tuple[float, float]:
    """Using Matlab Engine to do same thing as Pymap3d"""

    if has_mapping(eng):
        return eng.distance(lat1, lon1, lat2, lon2, eng.wgs84Ellipsoid(), nargout=2)
    else:
        return eng.matmap3d.vdist(lat1, lon1, lat2, lon2, nargout=2)


def stability(eng) -> bool:
    dlast, alast = nan, nan
    lon1, lon2 = 0.0, 1.0
    ok = True
    for i in range(20):
        lat1 = lat2 = 10.0 ** (-i)

        dist_m, az_deg = vdist(lat1, lon1, lat2, lon2)

        assert dist_m != dlast
        assert az_deg != alast
        dist_matlab, az_matlab = distance(eng, lat1, lon1, lat2, lon2)

        assert np.isreal(dist_m), f"Python vdist distance is not real for input latitude {lat1}"
        assert np.isreal(dist_matlab), f"Matlab distance is not real for input latitude {lat1}"

        if isclose(dist_matlab, dist_m) and isclose(az_matlab, az_deg):
            print(f"latitudes {lat1} {lat2}: {dist_m} meters {az_deg} deg azimuth")
        else:
            ok = False
            logging.error(
                f"""MISMATCH: latitude {lat1} {lat2}
    azimuth: Python: {az_matlab}  Matlab: {az_deg}
    distance: Python: {dist_m}  Matlab: {dist_matlab}
    """
            )

    return ok


eng = matlab_engine()

print("Mapping Toolbox:", has_mapping(eng))

if stability(eng):
    print("OK: vdist compare")
else:
    raise ValueError("FAIL: vdist compare")
