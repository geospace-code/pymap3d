import pytest
from pytest import approx
from datetime import datetime
import pymap3d as pm


def test_eci2ecef():
    pytest.importorskip("numpy")
    # this example from Matlab eci2ecef docs
    eci = [-2981784, 5207055, 3161595]
    utc = datetime(2019, 1, 4, 12)
    ecef = pm.eci2ecef(*eci, utc)

    rel = 0.0001
    assert ecef == approx([-5.7627e6, -1.6827e6, 3.1560e6], rel=rel)


def test_ecef2eci():
    pytest.importorskip("numpy")
    # this example from Matlab ecef2eci docs
    ecef = [-5762640, -1682738, 3156028]
    utc = datetime(2019, 1, 4, 12)
    eci = pm.ecef2eci(*ecef, utc)

    rel = 0.0001
    assert eci == approx([-2.9818e6, 5.2070e6, 3.1616e6], rel=rel)


def test_eci2geodetic():
    pytest.importorskip("numpy")

    eci = [-2981784, 5207055, 3161595]
    utc = datetime(2019, 1, 4, 12)
    lla = pm.eci2geodetic(*eci, utc)

    rel = 0.0001
    assert lla == approx([27.881, -163.722, 408850.65], rel=rel)


def test_geodetic2eci():
    pytest.importorskip("numpy")

    lla = [27.881, -163.722, 408850.65]
    utc = datetime(2019, 1, 4, 12)
    eci = pm.geodetic2eci(*lla, utc)

    rel = 0.0001
    assert eci == approx([-2981784, 5207055, 3161595], rel=rel)


def test_eci2aer():
    # test coords from Matlab eci2aer
    pytest.importorskip("numpy")
    t = datetime(1969, 7, 20, 21, 17, 40)

    eci = [-3.8454e8, -0.5099e8, -0.3255e8]
    lla = [28.4, -80.5, 2.7]

    aer = pm.eci2aer(*eci, *lla, t)
    rel = 0.0001
    assert aer == approx([162.55, 55.12, 384013940.9], rel=rel)


def test_aer2eci():
    # test coords from Matlab aer2eci
    pytest.importorskip("numpy")

    aer = [162.55, 55.12, 384013940.9]
    lla = [28.4, -80.5, 2.7]
    t = datetime(1969, 7, 20, 21, 17, 40)

    eci = pm.aer2eci(*aer, *lla, t)

    rel = 0.001
    assert eci == approx([-3.8454e8, -0.5099e8, -0.3255e8], rel=rel)

    with pytest.raises(ValueError):
        pm.aer2eci(aer[0], aer[1], -1, *lla, t)
