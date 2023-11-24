#!/usr/bin/env python3
"""
Compare with Matlab Mapping toolbox reckon()
"""

from __future__ import annotations
import argparse
import logging
from math import isclose

from matlab_engine import matlab_engine, has_mapping

from pymap3d.lox import loxodrome_direct


def reckon(eng, lat1: float, lon1: float, rng: float, az: float) -> tuple[float, float]:
    """Using Matlab Engine to do same thing as Pymap3d"""

    if has_mapping(eng):
        return eng.reckon("rh", lat1, lon1, rng, az, eng.wgs84Ellipsoid(), nargout=2)
    else:
        return eng.matmap3d.vreckon(lat1, lon1, rng, az, nargout=2)


def stability(eng) -> bool:
    clat, clon, rng = 35.0, 140.0, 50000.0  # arbitrary
    ok = True
    for i in range(20):
        for azi in (90 + 10.0 ** (-i), -90 + 10.0 ** (-i), 270 + 10.0 ** (-i), -270 + 10.0 ** (-i)):
            lat, lon = loxodrome_direct(clat, clon, rng, azi)

            lat_matlab, lon_matlab = reckon(eng, clat, clon, rng, azi)
            rstr = (
                f"azimuth: {azi} lat,lon: Python: {lat}, {lon}  Matlab: {lat_matlab}, {lon_matlab}"
            )
            if not (
                isclose(lat_matlab, lat, rel_tol=0.005) and isclose(lon_matlab, lon, rel_tol=0.001)
            ):
                ok = False
                logging.error(rstr)

    return ok


p = argparse.ArgumentParser(description="compare reckon loxodrome")
p.add_argument(
    "-f", "--force_matmap3d", help="use matmap3d instead of Matlab OEM Toolbox", action="store_true"
)
P = p.parse_args()

eng = matlab_engine()

print("Mapping Toolbox:", has_mapping(eng, P.force_matmap3d))

if stability(eng):
    print("OK: lox_stability: comparison")
else:
    raise ValueError("FAIL: lox_stability comparison")
