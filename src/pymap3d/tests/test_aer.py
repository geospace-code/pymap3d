from math import radians

import pymap3d as pm
import pytest
from pytest import approx

ELL = pm.Ellipsoid.from_name("wgs84")
A = ELL.semimajor_axis
B = ELL.semiminor_axis


@pytest.mark.parametrize(
    "aer,lla,xyz", [((33, 70, 1000), (42, -82, 200), (660930.2, -4701424.0, 4246579.6))]
)
def test_aer2ecef(aer, lla, xyz):
    # degrees
    xyz1 = pm.aer2ecef(*aer, *lla)
    assert xyz1 == approx(xyz)
    assert all(isinstance(n, float) for n in xyz1)
    # float includes np.float64 i.e. a scalar

    # radians
    raer = (radians(aer[0]), radians(aer[1]), aer[2])
    rlla = (radians(lla[0]), radians(lla[1]), lla[2])
    xyz1 = pm.aer2ecef(*raer, *rlla, deg=False)
    assert xyz1 == approx(xyz)
    assert all(isinstance(n, float) for n in xyz1)

    # bad input
    with pytest.raises(ValueError):
        pm.aer2ecef(aer[0], aer[1], -1, *lla)


@pytest.mark.parametrize(
    "xyz, lla, aer",
    [
        ((A - 1, 0, 0), (0, 0, 0), (0, -90, 1)),
        ((-A + 1, 0, 0), (0, 180, 0), (0, -90, 1)),
        ((0, A - 1, 0), (0, 90, 0), (0, -90, 1)),
        ((0, -A + 1, 0), (0, -90, 0), (0, -90, 1)),
        ((0, 0, B - 1), (90, 0, 0), (0, -90, 1)),
        ((0, 0, -B + 1), (-90, 0, 0), (0, -90, 1)),
        ((660930.19276, -4701424.22296, 4246579.60463), (42, -82, 200), (33, 70, 1000)),
    ],
)
def test_ecef2aer(xyz, lla, aer):
    # degrees
    aer1 = pm.ecef2aer(*xyz, *lla)
    assert aer1 == approx(aer)
    assert all(isinstance(n, float) for n in aer1)

    # radians
    rlla = (radians(lla[0]), radians(lla[1]), lla[2])
    raer = (radians(aer[0]), radians(aer[1]), aer[2])
    aer1 = pm.ecef2aer(*xyz, *rlla, deg=False)
    assert aer1 == approx(raer)
    assert all(isinstance(n, float) for n in aer1)


@pytest.mark.parametrize("aer,enu", [((33, 70, 1000), (186.2775, 286.8422, 939.6926))])
def test_aer_enu(aer, enu):
    # degrees
    enu1 = pm.aer2enu(*aer)
    assert enu1 == approx(enu)
    assert all(isinstance(n, float) for n in enu1)

    # radians
    raer = (radians(aer[0]), radians(aer[1]), aer[2])
    enu1 = pm.aer2enu(*raer, deg=False)
    assert enu1 == approx(enu)
    assert all(isinstance(n, float) for n in enu1)

    # bad input
    with pytest.raises(ValueError):
        pm.aer2enu(aer[0], aer[1], -1)

    # degrees
    aer1 = pm.enu2aer(*enu)
    assert aer1 == approx(aer)
    assert all(isinstance(n, float) for n in aer1)

    # radians
    aer1 = pm.enu2aer(*enu, deg=False)
    assert aer1 == approx(raer)
    assert all(isinstance(n, float) for n in aer1)


@pytest.mark.parametrize("aer,ned", [((33, 70, 1000), (286.8422, 186.2775, -939.6926))])
def test_aer_ned(aer, ned):
    ned1 = pm.aer2ned(*aer)
    assert ned1 == approx(ned)
    assert all(isinstance(n, float) for n in ned1)

    # bad value
    with pytest.raises(ValueError):
        pm.aer2ned(aer[0], aer[1], -1)

    aer1 = pm.ned2aer(*ned)
    assert aer1 == approx(aer)
    assert all(isinstance(n, float) for n in aer1)
