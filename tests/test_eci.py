import pytest
from pytest import approx
from datetime import datetime
import pymap3d as pm


def test_eci2ecef():

    pytest.importorskip("astropy")
    # this example from Matlab eci2ecef docs
    eci = [-2981784, 5207055, 3161595]
    utc = datetime(2019, 1, 4, 12)
    ecef = pm.eci2ecef(*eci, utc)
    assert ecef == approx([-5.7627e6, -1.6827e6, 3.1560e6], rel=0.01)


def test_ecef2eci():

    pytest.importorskip("astropy")
    # this example from Matlab ecef2eci docs
    ecef = [-5762640, -1682738, 3156028]
    utc = datetime(2019, 1, 4, 12)
    eci = pm.ecef2eci(*ecef, utc)
    print(eci)
    assert eci == approx([-2.9818e6, 5.2070e6, 3.1616e6], rel=0.01)


if __name__ == "__main__":
    pytest.main([__file__])
