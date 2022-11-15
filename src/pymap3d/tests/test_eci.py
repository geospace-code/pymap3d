from datetime import datetime

import pymap3d as pm
import pytest
from pytest import approx

try:
    import astropy
except ImportError:
    astropy = None


def test_eci2ecef() -> None:
    pytest.importorskip("numpy")
    # this example from Matlab eci2ecef docs
    eci = [-2981784.0, 5207055.0, 3161595.0]
    utc = datetime(2019, 1, 4, 12)
    ecef = pm.eci2ecef(*eci, utc)  # type: ignore[call-overload]

    rel = 0.025 if astropy is None else 0.0001

    assert ecef == approx([-5.7627e6, -1.6827e6, 3.1560e6], rel=rel)


def test_ecef2eci() -> None:
    pytest.importorskip("numpy")
    # this example from Matlab ecef2eci docs
    ecef = [-5762640, -1682738, 3156028]
    utc = datetime(2019, 1, 4, 12)
    eci = pm.ecef2eci(*ecef, utc)  # type: ignore[call-overload]

    rel = 0.01 if astropy is None else 0.0001

    assert eci == approx([-2981810.6, 5207039.5, 3161595.1], rel=rel)


def test_eci2geodetic() -> None:
    pytest.importorskip("numpy")

    eci = [-2981784, 5207055, 3161595]
    utc = datetime(2019, 1, 4, 12)
    lla = pm.eci2geodetic(*eci, utc)  # type: ignore[call-overload]

    rel = 0.01 if astropy is None else 0.0001

    assert lla == approx([27.880801, -163.722058, 408850.646], rel=rel)


def test_geodetic2eci() -> None:
    pytest.importorskip("numpy")

    lla = [27.880801, -163.722058, 408850.646]
    utc = datetime(2019, 1, 4, 12)
    eci = pm.geodetic2eci(*lla, utc)  # type: ignore[call-overload]

    rel = 0.01 if astropy is None else 0.0001

    assert eci == approx([-2981784, 5207055, 3161595], rel=rel)


def test_eci2aer() -> None:
    # test coords from Matlab eci2aer
    pytest.importorskip("numpy")
    t = datetime(1969, 7, 20, 21, 17, 40)

    eci = [-384535788.2, -51009524.3, -32566759.8]
    lla = [28.4, -80.5, 2.7]

    aer = pm.eci2aer(*eci, *lla, t)  # type: ignore[call-overload]

    rel = 0.01 if astropy is None else 0.0001

    assert aer == approx([162.55, 55.12, 384013940.9], rel=rel)


def test_aer2eci() -> None:
    # test coords from Matlab aer2eci
    pytest.importorskip("numpy")

    aer = [162.55, 55.12, 384013940.9]
    lla = [28.4, -80.5, 2.7]
    t = datetime(1969, 7, 20, 21, 17, 40)

    eci = pm.aer2eci(*aer, *lla, t)  # type: ignore[call-overload]

    rel = 0.1 if astropy is None else 0.001

    assert eci == approx([-384535788.2, -51009524.3, -32566759.8], rel=rel)

    with pytest.raises(ValueError):
        pm.aer2eci(aer[0], aer[1], -1, *lla, t)  # type: ignore[call-overload]
