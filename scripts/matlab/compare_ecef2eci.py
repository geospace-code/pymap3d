#!/usr/bin/env python3
"""
Compare ecef2eci() with Matlab Aerospace Toolbox
"""

from __future__ import annotations
import logging
from datetime import datetime

import numpy as np

import pymap3d

from matlab_engine import matlab_engine, has_aerospace, pydt2matdt, mat2np


def ecef2eci(eng, utc_m, ecef):
    if has_aerospace(eng):
        eci = eng.ecef2eci(utc_m, np.asarray(ecef), nargout=1)
        return mat2np(eci)

    return eng.matmap3d.ecef2eci(utc_m, *ecef, nargout=3)


def eci2ecef(eng, utc_m, eci):
    if has_aerospace(eng):
        ecef = eng.eci2ecef(utc_m, eci, nargout=1)
        return mat2np(ecef)

    return eng.matmap3d.eci2ecef(utc_m, *eci, nargout=3)


def compare_ecef2eci(eng, ecef, utc: datetime) -> bool:
    eci_py = pymap3d.ecef2eci(ecef[0], ecef[1], ecef[2], utc)

    eci_m = ecef2eci(eng, pydt2matdt(eng, utc), ecef)

    ok = bool(np.isclose(eci_py, eci_m, rtol=0.01).all())

    if not ok:
        logging.error(f"eci2ecef did not match Matlab  {eci_py}  !=  {eci_m}")

    return ok


def compare_eci2ecef(eng, eci, utc) -> bool:
    ecef_py = pymap3d.eci2ecef(eci[0], eci[1], eci[2], utc)

    ecef_m = eci2ecef(eng, pydt2matdt(eng, utc), eci)

    ok = bool(np.isclose(ecef_py, ecef_m, rtol=0.02).all())

    if not ok:
        logging.error(f"eci2ecef did not match Matlab  {ecef_py}  !=  {ecef_m}")

    return ok


eng = matlab_engine()

ecef = [-5762640.0, -1682738.0, 3156028.0]
eci = [-3009680.518620539, 5194367.153184303, 3156028.0]
utc_py = datetime(2019, 1, 4, 12)

print("Aerospace Toolbox:", has_aerospace(eng))

if ecef2eci_ok := compare_ecef2eci(eng, ecef, utc_py):
    print("OK: PyMap3d ecef2eci vs. Matlab ecef2eci")

if eci2ecef_ok := compare_eci2ecef(eng, eci, utc_py):
    print("OK: PyMap3d eci2ecef vs. Matlab eci2ecef")

if not ecef2eci_ok and eci2ecef_ok:
    raise ValueError("compare_ecef2eci: Matlab compare mismatch")
