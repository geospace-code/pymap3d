#!/usr/bin/env python3
"""
Compare with Matlab Mapping Toolbox reckon()
"""

from __future__ import annotations

import logging
from math import isclose, nan
import numpy as np

from matlab_engine import matlab_engine, has_mapping

import pymap3d.vincenty


def reckon(eng, lat1: float, lon1: float, srng: float, az: float) -> tuple[float, float]:
    """Using Matlab Engine to do same thing as Pymap3d"""

    if has_mapping(eng):
        return eng.reckon("gc", lat1, lon1, srng, az, eng.wgs84Ellipsoid(), nargout=2)
    else:
        return eng.matmap3d.vreckon(lat1, lon1, srng, az, nargout=2)


def stability(eng) -> bool:
    dlast, alast = nan, nan
    lon1, lon2 = 0.0, 1.0
    ok = True
    for i in range(20):
        lat1 = lat2 = 10.0 ** (-i)

        dist_m, az_deg = pymap3d.vincenty.vreckon(lat1, lon1, lat2, lon2)

        assert dist_m != dlast
        assert az_deg != alast
        dist_matlab, az_matlab = reckon(eng, lat1, lon1, lat2, lon2)

        assert np.isreal(dist_m), f"Python reckon is not real for input latitude {lat1}"
        assert np.isreal(dist_matlab), f"Matlab reckon is not real for input latitude {lat1}"

        if isclose(dist_matlab, dist_m) and isclose(az_matlab, az_deg, rel_tol=0.005):
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
    print("OK: vreckon compare")
else:
    raise ValueError("FAIL: vreckon compare")
