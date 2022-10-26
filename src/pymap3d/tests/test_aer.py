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
    x, y, z = pm.aer2ecef(*aer, *lla)
    assert x == approx(xyz[0])
    assert y == approx(xyz[1])
    assert z == approx(xyz[2])
    assert isinstance(x, float)
    assert isinstance(y, float)
    assert isinstance(z, float)

    raer = (radians(aer[0]), radians(aer[1]), aer[2])
    rlla = (radians(lla[0]), radians(lla[1]), lla[2])
    assert pm.aer2ecef(*raer, *rlla, deg=False) == approx(xyz)

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
    assert pm.ecef2aer(*xyz, *lla) == approx(aer)

    rlla = (radians(lla[0]), radians(lla[1]), lla[2])
    raer = (radians(aer[0]), radians(aer[1]), aer[2])
    assert pm.ecef2aer(*xyz, *rlla, deg=False) == approx(raer)


@pytest.mark.parametrize("aer,enu", [((33, 70, 1000), (186.2775, 286.8422, 939.6926))])
def test_aer_enu(aer, enu):
    e, n, u = pm.aer2enu(*aer)
    assert e == approx(enu[0])
    assert n == approx(enu[1])
    assert u == approx(enu[2])
    assert isinstance(e, float)
    assert isinstance(n, float)
    assert isinstance(u, float)

    raer = (radians(aer[0]), radians(aer[1]), aer[2])
    assert pm.aer2enu(*raer, deg=False) == approx(enu)

    with pytest.raises(ValueError):
        pm.aer2enu(aer[0], aer[1], -1)

    a, e, r = pm.enu2aer(*enu)
    assert a == approx(aer[0])
    assert e == approx(aer[1])
    assert r == approx(aer[2])
    assert isinstance(a, float)
    assert isinstance(e, float)
    assert isinstance(r, float)

    assert pm.enu2aer(*enu, deg=False) == approx(raer)


@pytest.mark.parametrize("aer,ned", [((33, 70, 1000), (286.8422, 186.2775, -939.6926))])
def test_aer_ned(aer, ned):
    assert pm.aer2ned(*aer) == approx(ned)

    with pytest.raises(ValueError):
        pm.aer2ned(aer[0], aer[1], -1)

    assert pm.ned2aer(*ned) == approx(aer)
