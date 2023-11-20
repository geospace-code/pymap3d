#!/usr/bin/env python3
from __future__ import annotations

from datetime import datetime

from pytest import approx

from pymap3d.eci import ecef2eci, eci2ecef

from matlab_aerospace import matlab_aerospace
from matlab_engine import matlab_engine


def test_ecef_eci():
    ecef = [-5762640.0, -1682738.0, 3156028.0]
    utc = datetime(2019, 1, 4, 12)

    eci = ecef2eci(*ecef, utc)

    eng = matlab_engine()

    utc_matlab = eng.datetime(utc.year, utc.month, utc.day, utc.hour, utc.minute, utc.second)
    if matlab_aerospace(eng):
        eci_matlab = eng.ecef2eci(utc_matlab, *ecef, nargout=3)
    else:
        eci_matlab = eng.matmap3d.ecef2eci(utc_matlab, *ecef, nargout=3)

    assert eci == approx(eci_matlab, rel=0.01)

    print("OK: PyMap3d ecef2eci vs. Matlab ecef2eci")

    ecef = eci2ecef(*eci_matlab, utc)

    if matlab_aerospace(eng):
        ecef_matlab = eng.eci2ecef(utc_matlab, *eci_matlab, nargout=3)
    else:
        ecef_matlab = eng.matmap3d.eci2ecef(utc_matlab, *eci_matlab, nargout=3)

    assert ecef == approx(ecef_matlab, rel=0.02)

    print("OK: PyMap3d eci2ecef vs. Matlab eci2ecef")


if __name__ == "__main__":
    test_ecef_eci()
