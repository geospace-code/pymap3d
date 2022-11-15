#!/usr/bin/env python3
from __future__ import annotations

from math import isclose
from pathlib import Path
from datetime import datetime

import matlab.engine
from pymap3d.eci import ecef2eci

cwd = Path(__file__).parent
eng = None  # don't start Matlab engine over and over when script is interactive

if eng is None:
    eng = matlab.engine.start_matlab("-nojvm")
    eng.addpath(eng.genpath(str(cwd)), nargout=0)

if eng.has_aerospace_toolbox():
    has_aero = True
else:
    has_aero = False
    d = cwd.parents[1] / "matmap3d"
    if d.is_dir():
        eng.addpath(str(d), nargout=0)
    else:
        raise EnvironmentError(f"Matlab {eng.version()} does not have Aerospace Toolbox")


ecef = [-5762640, -1682738, 3156028]
utc = datetime(2019, 1, 4, 12)

eci = ecef2eci(*ecef, utc)

utc_matlab = eng.datetime(utc.year, utc.month, utc.day, utc.hour, utc.minute, utc.second)
if has_aero:
    eci_matlab = eng.ecef2eci(utc_matlab, *ecef, nargout=3)  # type: ignore
else:
    eci_matlab = eng.matmap3d.ecef2eci(utc_matlab, *ecef, nargout=3)  # type: ignore

assert isclose(eci[0], eci_matlab[0], rel_tol=0.01)
assert isclose(eci[1], eci_matlab[1], rel_tol=0.01)
assert isclose(eci[2], eci_matlab[2], rel_tol=0.01)

print("OK: PyMap3d ecef2eci vs. Matlab ecef2eci")
