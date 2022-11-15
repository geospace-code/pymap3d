from datetime import datetime

import pymap3d as pm
import pytest
from pytest import approx

try:
    import astropy
except ImportError:
    astropy = None


def test_eci2ecef():
    pytest.importorskip("numpy")
    # this example from Matlab eci2ecef docs
    eci = [-2981784, 5207055, 3161595]
    utc = datetime(2019, 1, 4, 12)
    ecef = pm.eci2ecef(*eci, utc)

    rel = 0.025 if astropy is None else 0.0001

    assert ecef == approx([-5.7627e6, -1.6827e6, 3.1560e6], rel=rel)


def test_ecef2eci():
    pytest.importorskip("numpy")
    # this example from Matlab ecef2eci docs
    ecef = [-5762640, -1682738, 3156028]
    utc = datetime(2019, 1, 4, 12)
    eci = pm.ecef2eci(*ecef, utc)

    rel = 0.01 if astropy is None else 0.0001

    assert eci == approx([-2981810.6, 5207039.5, 3161595.1], rel=rel)


def test_eci2geodetic():
    pytest.importorskip("numpy")

    eci = [-2981784, 5207055, 3161595]
    utc = datetime(2019, 1, 4, 12)
    lla = pm.eci2geodetic(*eci, utc)

    rel = 0.01 if astropy is None else 0.0001

    assert lla == approx([27.880801, -163.722058, 408850.646], rel=rel)


def test_geodetic2eci():
    pytest.importorskip("numpy")

    lla = [27.880801, -163.722058, 408850.646]
    utc = datetime(2019, 1, 4, 12)
    eci = pm.geodetic2eci(*lla, utc)

    rel = 0.01 if astropy is None else 0.0001

    assert eci == approx([-2981784, 5207055, 3161595], rel=rel)


def test_eci_aer():
    # test coords from Matlab eci2aer
    pytest.importorskip("numpy")
    t = datetime(2022, 1, 2, 3, 4, 5)

    eci = [4500000, -45000000, 3000000]
    lla = [28, -80, 100]

    aer = pm.eci2aer(*eci, *lla, t)

    rel = 0.01 if astropy is None else 0.0001

    assert aer == approx([314.9945, -53.0089, 5.026e7], rel=rel)

    eci2 = pm.aer2eci(*aer, *lla, t)

    rel = 0.1 if astropy is None else 0.001

    assert eci2 == approx(eci, rel=rel)

    with pytest.raises(ValueError):
        pm.aer2eci(aer[0], aer[1], -1, *lla, t)
