"""
Compare with Matlab Mapping Toolbox distance()
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


from pymap3d.vincenty import vdist


def distance(eng, matmap3d: bool, lat1, lon1, lat2, lon2) -> tuple[float, float]:
    """Using Matlab Engine to do same thing as Pymap3d"""

    if matmap3d:
        return eng.matmap3d.vdist(lat1, lon1, lat2, lon2, nargout=2)

    return eng.distance(lat1, lon1, lat2, lon2, eng.wgs84Ellipsoid(), nargout=2)


@pytest.mark.parametrize("matmap3d", [False, True])
def test_matlab_stability(matmap3d):
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

        dist_m, az_deg = vdist(lat1, lon1, lat2, lon2)

        assert dist_m != dlast
        assert az_deg != alast
        dist_matlab, az_matlab = distance(eng, matmap3d, lat1, lon1, lat2, lon2)

        assert dist_m == approx(dist_matlab)
        assert az_deg == approx(az_matlab)
