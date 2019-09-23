#!/usr/bin/env python
import pytest
from pytest import approx
import pymap3d as pm
from math import nan

lla0 = (42, -82, 200)
az = [0.0, 10.0, 125.0]


def test_losint():

    lat, lon, sr = pm.lookAtSpheroid(*lla0, az[1], tilt=0.0)

    assert (lat, lon, sr) == approx(lla0)

    with pytest.raises(ValueError):
        pm.lookAtSpheroid(lla0[0], lla0[1], -1, az, 0)


def test_array():
    np = pytest.importorskip("numpy")
    tilt = [30.0, 45.0, 90.0]
    lat, lon, sr = pm.lookAtSpheroid(*lla0, az, tilt)

    truth = np.array([[42.00103959, lla0[1], 230.9413173], [42.00177328, -81.9995808, 282.84715651], [nan, nan, nan]])

    assert np.column_stack((lat, lon, sr)) == approx(truth, nan_ok=True)
