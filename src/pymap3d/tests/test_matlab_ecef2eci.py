"""
Compare ecef2eci() with Matlab Aerospace Toolbox
"""

from __future__ import annotations
from datetime import datetime

import pytest
from pytest import approx

import pymap3d

try:
    import numpy as np
    from .matlab_engine import matlab_engine, has_aerospace, has_matmap3d, pydt2matdt
except ImportError:
    pytest.skip("Matlab Engine not found", allow_module_level=True)


def ecef2eci(eng, matmap3d: bool, utc_m, ecef):
    if matmap3d:
        return eng.matmap3d.ecef2eci(utc_m, *ecef, nargout=3)

    return np.array(eng.ecef2eci(utc_m, np.asarray(ecef), nargout=1)).squeeze()


def eci2ecef(eng, matmap3d: bool, utc_m, eci):
    if matmap3d:
        return eng.matmap3d.eci2ecef(utc_m, *eci, nargout=3)

    return np.array(eng.eci2ecef(utc_m, np.asarray(eci), nargout=1)).squeeze()


@pytest.mark.parametrize("matmap3d", [False, True])
def test_compare_ecef2eci(matmap3d):
    eng = matlab_engine()

    if matmap3d:
        if not has_matmap3d(eng):
            pytest.skip("Matmap3d not found")
    else:
        if not has_aerospace(eng):
            pytest.skip("Aerospace Toolbox not found")

    ecef = [-5762640.0, -1682738.0, 3156028.0]
    utc = datetime(2019, 1, 4, 12)
    rtol = 0.01

    eci_py = pymap3d.ecef2eci(ecef[0], ecef[1], ecef[2], utc)

    eci_m = ecef2eci(eng, matmap3d, pydt2matdt(eng, utc), ecef)

    assert eci_py == approx(eci_m, rel=rtol)


@pytest.mark.parametrize("matmap3d", [False, True])
def test_compare_eci2ecef(matmap3d):
    eng = matlab_engine()

    if matmap3d:
        if not has_matmap3d(eng):
            pytest.skip("Matmap3d not found")
    else:
        if not has_aerospace(eng):
            pytest.skip("Aerospace Toolbox not found")

    eci = [-3009680.518620539, 5194367.153184303, 3156028.0]
    utc = datetime(2019, 1, 4, 12)
    rtol = 0.02

    ecef_py = pymap3d.eci2ecef(eci[0], eci[1], eci[2], utc)

    ecef_m = eci2ecef(eng, matmap3d, pydt2matdt(eng, utc), eci)

    assert ecef_py == approx(ecef_m, rel=rtol)
