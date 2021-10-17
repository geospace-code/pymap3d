#!/usr/bin/env python3
from __future__ import annotations
import logging
from pathlib import Path
from math import isclose, nan

from pymap3d.lox import loxodrome_direct

import matlab.engine

cwd = Path(__file__).parent
eng = None  # don't start Matlab engine over and over when script is interactive

if eng is None:
    eng = matlab.engine.start_matlab("-nojvm")
    eng.addpath(eng.genpath(str(cwd)), nargout=0)

if not eng.has_map_toolbox():
    raise EnvironmentError("Matlab does not have Mapping Toolbox")


def matlab_func(lat1: float, lon1: float, rng: float, az: float) -> tuple[float, float]:
    """Using Matlab Engine to do same thing as Pymap3d"""
    return eng.reckon("rh", lat1, lon1, rng, az, eng.wgs84Ellipsoid(), nargout=2)  # type: ignore


lat_last, lon_last = nan, nan
clat, clon, rng = 40.0, -80.0, 10.0  # arbitrary

for i in range(20):
    azi = 90.0 + 10.0 ** (-i)

    lat, lon = loxodrome_direct(clat, clon, rng, azi)

    assert lat != lat_last
    assert lon != lon_last

    lat_matlab, lon_matlab = matlab_func(clat, clon, rng, azi)
    rstr = f"azimuth: {azi} lat,lon: Python: {lat}, {lon}  Matlab: {lat_matlab}, {lon_matlab}"
    if not (isclose(lat_matlab, lat, rel_tol=0.005) and isclose(lon_matlab, lon, rel_tol=0.001)):
        logging.error(rstr)
