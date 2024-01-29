"""
Compare with Matlab Mapping Toolbox reckon()
"""

from __future__ import annotations
from math import nan

import pytest
from pytest import approx

try:
    from .matlab_engine import matlab_engine, has_mapping, has_matmap3d
except ImportError:
    pytest.skip("Matlab Engine not found", allow_module_level=True)
except RuntimeError:
    pytest.skip("Matlab Engine configuration error", allow_module_level=True)


import pymap3d.vincenty


def reckon(
    eng, matmap3d: bool, lat1: float, lon1: float, srng: float, az: float
) -> tuple[float, float]:
    """Using Matlab Engine to do same thing as Pymap3d"""

    if matmap3d:
        return eng.matmap3d.vreckon(lat1, lon1, srng, az, nargout=2)

    return eng.reckon("gc", lat1, lon1, srng, az, eng.wgs84Ellipsoid(), nargout=2)


@pytest.mark.parametrize("matmap3d", [False, True])
def test_reckon_stability(matmap3d):
    eng = matlab_engine()

    if matmap3d:
        if not has_matmap3d(eng):
            pytest.skip("Matmap3d not found")
    else:
        if not has_mapping(eng):
            pytest.skip("Matlab Toolbox not found")

    dlast, alast = nan, nan
    lon1, lon2 = 0.0, 1.0
    for i in range(20):
        lat1 = lat2 = 10.0 ** (-i)

        dist_m, az_deg = pymap3d.vincenty.vreckon(lat1, lon1, lat2, lon2)

        assert dist_m != dlast
        assert az_deg != alast
        dist_matlab, az_matlab = reckon(eng, matmap3d, lat1, lon1, lat2, lon2)

        assert dist_m == approx(dist_matlab)

        assert az_deg == approx(az_matlab, rel=0.005)


@pytest.mark.parametrize("matmap3d", [False, True])
def test_reckon_unit(matmap3d):
    """
    Test various extrema and other values of interest
    """

    eng = matlab_engine()

    if matmap3d:
        if not has_matmap3d(eng):
            pytest.skip("Matmap3d not found")
    else:
        if not has_mapping(eng):
            pytest.skip("Matlab Toolbox not found")

    latlon88 = 52.22610277777778, -1.2696583333333333
    srng88 = 839.63
    az88 = 63.02

    # issue 88
    lat_p, lon_p = pymap3d.vincenty.vreckon(*latlon88, srng88, az88)

    lat_m, lon_m = reckon(eng, matmap3d, *latlon88, srng88, az88)

    assert lat_p == approx(lat_m)
    assert lon_p == approx(lon_m)
