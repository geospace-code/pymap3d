#!/usr/bin/env python3
from __future__ import annotations

import logging
from math import isclose
from pathlib import Path

from matlab_mapping import matlab_mapping

import matlab.engine
from pymap3d.lox import loxodrome_direct


cwd = Path(__file__).parent
eng = matlab.engine.start_matlab("-nojvm")
eng.addpath(eng.genpath(str(cwd)), nargout=0)

has_map = matlab_mapping(eng)


def reckon(lat1: float, lon1: float, rng: float, az: float) -> tuple[float, float]:
    """Using Matlab Engine to do same thing as Pymap3d"""
    if has_map:
        return eng.reckon("rh", lat1, lon1, rng, az, eng.wgs84Ellipsoid(), nargout=2)
    else:
        return eng.matmap3d.vreckon(lat1, lon1, rng, az, nargout=2)


clat, clon, rng = 35.0, 140.0, 50000.0  # arbitrary

Nerr = 0
for i in range(20):
    for azi in (90 + 10.0 ** (-i), -90 + 10.0 ** (-i), 270 + 10.0 ** (-i), -270 + 10.0 ** (-i)):
        lat, lon = loxodrome_direct(clat, clon, rng, azi)

        lat_matlab, lon_matlab = reckon(clat, clon, rng, azi)
        rstr = f"azimuth: {azi} lat,lon: Python: {lat}, {lon}  Matlab: {lat_matlab}, {lon_matlab}"
        if not (
            isclose(lat_matlab, lat, rel_tol=0.005) and isclose(lon_matlab, lon, rel_tol=0.001)
        ):
            Nerr += 1
            logging.error(rstr)

if Nerr == 0:
    print("lox_stability: comparison OK")
