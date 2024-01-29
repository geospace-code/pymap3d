#!/usr/bin/env python3
"""
Compare with Matlab Mapping toolbox reckon()
"""

from __future__ import annotations

import pytest
from pytest import approx

try:
    import numpy as np
    from .matlab_engine import matlab_engine, has_mapping
except ImportError:
    pytest.skip("Matlab Engine not found", allow_module_level=True)
except RuntimeError:
    pytest.skip("Matlab Engine configuration error", allow_module_level=True)


import pymap3d.vincenty


def track2(eng, lat1: float, lon1: float, lat2: float, lon2: float, npts: int, deg: bool) -> tuple:
    """Using Matlab Engine to do same thing as Pymap3d"""
    d = "degrees" if deg else "radians"

    lats, lons = eng.track2(
        "gc", lat1, lon1, lat2, lon2, eng.wgs84Ellipsoid(), d, float(npts), nargout=2
    )
    return np.array(lats).squeeze(), np.array(lons).squeeze()


@pytest.mark.parametrize("deg", [True, False])
def test_track2_compare(deg):
    lat1, lon1 = 0.0, 80.0
    lat2, lon2 = 0.0, 81.0
    if not deg:
        lat1, lon1, lat2, lon2 = np.radians((lat1, lon1, lat2, lon2))

    eng = matlab_engine()

    if not has_mapping(eng):
        pytest.skip("Matlab Toolbox not found")

    lats, lons = pymap3d.vincenty.track2(lat1, lon1, lat2, lon2, npts=4, deg=deg)

    lats_m, lons_m = track2(eng, lat1, lon1, lat2, lon2, npts=4, deg=deg)

    assert lats == approx(lats_m)
    assert lons == approx(lons_m)
