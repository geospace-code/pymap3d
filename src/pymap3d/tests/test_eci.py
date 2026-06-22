import datetime

import pymap3d as pm
import pytest
from pytest import approx

try:
    import astropy
except ImportError:
    astropy = None

ECI = (-2981784.0, 5207055.0, 3161595.0)
ECEF = [-5762640.0, -1682738.0, 3156028.0]
UTC = datetime.datetime(2019, 1, 4, 12, tzinfo=datetime.timezone.utc)


def test_eci2ecef():
    pytest.importorskip("numpy")
    # this example from Matlab eci2ecef docs
    ecef = pm.eci2ecef(*ECI, UTC)

    assert isinstance(ecef[0], float)
    assert isinstance(ecef[1], float)
    assert isinstance(ecef[2], float)


def test_eci2ecef_numpy():
    ecef = pm.eci2ecef(*ECI, UTC, force_non_astropy=True)

    assert ecef == approx(ECEF, abs=60)
    assert isinstance(ecef[0], float)
    assert isinstance(ecef[1], float)
    assert isinstance(ecef[2], float)


def test_eci2ecef_astropy():
    pytest.importorskip("astropy")

    ecef = pm.eci2ecef(*ECI, UTC)

    rel = 0.0001

    assert ecef == approx(ECEF, rel=rel)
    assert isinstance(ecef[0], float)
    assert isinstance(ecef[1], float)
    assert isinstance(ecef[2], float)


@pytest.mark.parametrize("force_non_astropy", [True, False])
def test_ecef2eci(force_non_astropy):
    pytest.importorskip("numpy")
    # this example from Matlab ecef2eci docs
    eci = pm.ecef2eci(*ECEF, UTC, force_non_astropy=force_non_astropy)

    assert isinstance(eci[0], float)
    assert isinstance(eci[1], float)
    assert isinstance(eci[2], float)


def test_ecef2eci_numpy():
    eci = pm.eci.ecef2eci_numpy(*ECEF, UTC)

    assert eci == approx(ECI, abs=60)
    assert isinstance(eci[0], float)
    assert isinstance(eci[1], float)
    assert isinstance(eci[2], float)


def test_ecef2eci_astropy():
    pytest.importorskip("astropy")

    eci = pm.eci.ecef2eci_astropy(*ECEF, UTC)

    rel = 0.0001

    assert eci == approx(ECI, rel=rel)
    assert isinstance(eci[0], float)
    assert isinstance(eci[1], float)
    assert isinstance(eci[2], float)


def test_eci2geodetic():
    lla = pm.eci2geodetic(*ECI, UTC)

    rel = 2e-6 if astropy is None else 0.0001

    assert lla == approx([27.880801, -163.722058, 408850.646], rel=rel)


def test_geodetic2eci():
    lla = [27.880801, -163.722058, 408850.646]

    eci = pm.geodetic2eci(*lla, UTC)

    rel = 2e-5 if astropy is None else 0.0001

    assert eci == approx([-2981784, 5207055, 3161595], rel=rel)


def test_eci_ecef_eop_roundtrip():
    eop = {"delta_ut1": -0.139198, "xp": 0.094, "yp": 0.324}

    ecef = pm.eci2ecef(*ECI, UTC, force_non_astropy=True, **eop)
    eci = pm.ecef2eci(*ecef, UTC, force_non_astropy=True, **eop)

    assert eci == approx(ECI, abs=1e-6)


def test_eci_aer():
    # test coords from Matlab eci2aer
    pytest.importorskip("numpy")
    t = datetime.datetime(2022, 1, 2, 3, 4, 5, tzinfo=datetime.timezone.utc)

    eci = [4500000, -45000000, 3000000]
    lla = [28, -80, 100]

    aer = pm.eci2aer(*eci, *lla, t)

    rel = 0.01 if astropy is None else 0.0001

    assert aer == approx([314.9945, -53.0089, 5.026e7], rel=rel)

    eci2 = pm.aer2eci(*aer, *lla, t)

    rel = 0.1 if astropy is None else 0.001

    assert eci2 == approx(eci, rel=rel)
