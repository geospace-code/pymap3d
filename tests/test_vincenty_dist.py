#!/usr/bin/env python3
import pytest
from pytest import approx
import pymap3d.vincenty as vincenty
import numpy as np


@pytest.mark.parametrize(
    "lat,lon,lat1,lon1,srange,az,backaz",
    [
        (0, 0, 0, 0, 0, 0, 0),
        (0, 0, 0, 90, 1.001875e7, 90, 270),
        (0, 0, 0, -90, 1.001875e7, 270, 90),
        (0, 0, 0, 180, 2.00375e7, 90, 270),
        (0, 0, 0, -180, 2.00375e7, 90, 270),
        (0, 0, 0, 4, 445277.96, 90, 270),
        (0, 0, 0, 5, 556597.45, 90, 270),
        (0, 0, 0, 6, 667916.94, 90, 270),
        (0, 0, 0, -6, 667916.94, 270, 90),
        (0, 0, 0, 7, 779236.44, 90, 270),
        (1e-16, 1e-16, 1e-16, 1, 111319.49, 90, 270),
        (90, 0, 0, 0, 1.00019657e7, 180, 0),
        (90, 0, -90, 0, 2.000393145e7, 180, 0),
    ],
)
def test_unit(lat, lon, lat1, lon1, srange, az, backaz):
    dist, az1, backaz1 = vincenty.vdist(lat, lon, lat1, lon1)
    assert dist == approx(srange, rel=0.005)
    assert az1 == approx(az)
    assert backaz1 == approx(backaz)


def test_vector():
    asr, aaz, aa21 = vincenty.vdist(10, 20, [10.02137267, 10.01917819], [20.0168471, 20.0193493])

    assert np.all(3e3 == approx(asr))  # for older pytest
    assert aaz == approx([38, 45])


@pytest.mark.parametrize("lat,lon,slantrange,az", [(10, 20, 3e3, 38), (0, 0, 0, 0)])
def test_identity(lat, lon, slantrange, az):
    lat1, lon1, a21 = vincenty.vreckon(lat, lon, slantrange, az)

    dist, az1, backaz1 = vincenty.vdist(lat, lon, lat1, lon1)

    assert dist == approx(slantrange)
    assert az1 == approx(az)


if __name__ == "__main__":
    pytest.main(["-v", __file__])
