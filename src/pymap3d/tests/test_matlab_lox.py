#!/usr/bin/env python3
"""
Compare with Matlab Mapping toolbox reckon()
"""

from __future__ import annotations

import pytest
from pytest import approx

try:
    from .matlab_engine import matlab_engine, has_matmap3d, has_mapping
except ImportError:
    pytest.skip("Matlab Engine not found", allow_module_level=True)


from pymap3d.lox import loxodrome_direct


def reckon(
    eng, matmap3d: bool, lat1: float, lon1: float, rng: float, az: float
) -> tuple[float, float]:
    """Using Matlab Engine to do same thing as Pymap3d"""

    if matmap3d:
        return eng.matmap3d.vreckon(lat1, lon1, rng, az, nargout=2)

    return eng.reckon("rh", lat1, lon1, rng, az, eng.wgs84Ellipsoid(), nargout=2)


@pytest.mark.parametrize("matmap3d", [False, True])
def test_lox_stability(matmap3d):
    eng = matlab_engine()

    if matmap3d:
        if not has_matmap3d(eng):
            pytest.skip("Matmap3d not found")
    else:
        if not has_mapping(eng):
            pytest.skip("Matlab Toolbox not found")

    clat, clon, rng = 35.0, 140.0, 50000.0  # arbitrary

    for i in range(20):
        for azi in (90 + 10.0 ** (-i), -90 + 10.0 ** (-i), 270 + 10.0 ** (-i), -270 + 10.0 ** (-i)):
            lat, lon = loxodrome_direct(clat, clon, rng, azi)

            lat_matlab, lon_matlab = reckon(eng, matmap3d, clat, clon, rng, azi)

            assert lat == approx(lat_matlab, rel=0.005)
            assert lon == approx(lon_matlab, rel=0.001)
