#!/usr/bin/env python3
from math import radians
import pytest
from pytest import approx

import pymap3d as pm

ELL = pm.Ellipsoid()
A = ELL.semimajor_axis
B = ELL.semiminor_axis


@pytest.mark.parametrize("xyz", [(0, A, 50), ([0], [A], [50])], ids=("scalar", "list"))
def test_scalar_enu(xyz):
    """
    verify we can handle the wide variety of input data type users might use
    """
    if isinstance(xyz[0], list):
        pytest.importorskip("numpy")

    enu = pm.ecef2enu(*xyz, 0, 90, -100)

    assert pm.enu2ecef(*enu, 0, 90, -100) == approx([0, A, 50])


def test_3d_enu():
    np = pytest.importorskip("numpy")
    xyz = (np.atleast_3d(0), np.atleast_3d(A), np.atleast_3d(50))

    enu = pm.ecef2enu(*xyz, 0, 90, -100)
    assert pm.enu2ecef(*enu, 0, 90, -100) == approx([0, A, 50])


@pytest.mark.parametrize("enu,lla,xyz", [((0, 0, 0), (0, 0, 0), (A, 0, 0)), ((0, 0, 1000), (0, 0, 0), (A + 1000, 0, 0))])
def test_enu_ecef(enu, lla, xyz):
    x, y, z = pm.enu2ecef(*enu, *lla)
    assert x == approx(xyz[0])
    assert y == approx(xyz[1])
    assert z == approx(xyz[2])
    assert isinstance(x, float)
    assert isinstance(y, float)
    assert isinstance(z, float)

    rlla = (radians(lla[0]), radians(lla[1]), lla[2])
    assert pm.enu2ecef(*enu, *rlla, deg=False) == approx(xyz)

    e, n, u = pm.ecef2enu(*xyz, *lla)
    assert e == approx(enu[0])
    assert n == approx(enu[1])
    assert u == approx(enu[2])
    assert isinstance(e, float)
    assert isinstance(n, float)
    assert isinstance(u, float)

    e, n, u = pm.ecef2enu(*xyz, *rlla, deg=False)
    assert e == approx(enu[0])
    assert n == approx(enu[1])
    assert u == approx(enu[2])


if __name__ == "__main__":
    pytest.main([__file__])
