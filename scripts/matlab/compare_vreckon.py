#!/usr/bin/env python3
"""
Compare with Matlab Mapping Toolbox reckon()
"""

from __future__ import annotations
import argparse
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


def unit(eng) -> bool:
    """
    Test various extrema and other values of interest
    """

    latlon88 = 52.22610277777778, -1.2696583333333333
    srng88 = 839.63
    az88 = 63.02

    # issue 88
    lat_p, lon_p = pymap3d.vincenty.vreckon(*latlon88, srng88, az88)

    lat_m, lon_m = reckon(eng, *latlon88, srng88, az88)

    ok = isclose(lat_p, lat_m) and isclose(lon_p, lon_m)

    if not ok:
        logging.error(
            f"""MISMATCH: lat/lon {latlon88[0]} {latlon88[1]}
azimuth: Python: {az88}  Matlab: {az88}
distance: Python: {srng88}  Matlab: {srng88}
"""
        )

    return ok


p = argparse.ArgumentParser(description="compare vreckon")
p.add_argument(
    "-f", "--force_matmap3d", help="use matmap3d instead of Matlab OEM Toolbox", action="store_true"
)
P = p.parse_args()

eng = matlab_engine()

print("Mapping Toolbox:", has_mapping(eng, P.force_matmap3d))

if unit(eng) and stability(eng):
    print("OK: vreckon compare")
else:
    raise ValueError("FAIL: vreckon compare")
