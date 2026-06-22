import pymap3d.vincenty as vincenty

import pytest
from pytest import approx

from math import radians


@pytest.mark.parametrize("deg", [True, False])
def test_track2_unit(deg):
    np = pytest.importorskip("numpy")

    lat1, lon1 = 0.0, 80.0
    lat2, lon2 = 0.0, 81.0
    lat0 = [0.0, 0.0, 0.0, 0.0]
    lon0 = [80.0, 80.33333, 80.66666, 81.0]
    if not deg:
        lat1 = radians(lat1)
        lon1 = radians(lon1)
        lat2 = radians(lat2)
        lon2 = radians(lon2)
        lon0 = list(map(radians, lon0))

    lats, lons = vincenty.track2(lat1, lon1, lat2, lon2, npts=4, deg=deg)

    assert lats == approx(lat0)
    assert lons == approx(lon0)
