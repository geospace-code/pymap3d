#!/usr/bin/env python3
from math import radians
import pytest
from pytest import approx
import pymap3d as pm

ELL = pm.Ellipsoid()
A = ELL.semimajor_axis
B = ELL.semiminor_axis


@pytest.mark.parametrize("aer,lla,xyz", [((33, 70, 1000), (42, -82, 200), (660930.2, -4701424.0, 4246579.6))])
def test_aer2ecef(aer, lla, xyz):
    x, y, z = pm.aer2ecef(*aer, *lla)
    assert x == approx(xyz[0])
    assert y == approx(xyz[1])
    assert z == approx(xyz[2])

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
    assert pm.aer2enu(*aer) == approx(enu)

    raer = (radians(aer[0]), radians(aer[1]), aer[2])
    assert pm.aer2enu(*raer, deg=False) == approx(enu)

    with pytest.raises(ValueError):
        pm.aer2enu(aer[0], aer[1], -1)

    assert pm.enu2aer(*enu) == approx(aer)
    assert pm.enu2aer(*enu, deg=False) == approx(raer)


@pytest.mark.parametrize("aer,ned", [((33, 70, 1000), (286.8422, 186.2775, -939.6926))])
def test_ned(aer, ned):
    assert pm.aer2ned(*aer) == approx(ned)

    with pytest.raises(ValueError):
        pm.aer2ned(aer[0], aer[1], -1)

    assert pm.ned2aer(*ned) == approx(aer)


if __name__ == "__main__":
    pytest.main(["-v", __file__])
