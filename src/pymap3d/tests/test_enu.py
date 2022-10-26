from math import radians

import pymap3d as pm
import pytest
from pytest import approx

ELL = pm.Ellipsoid.from_name("wgs84")
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
    assert pm.enu2ecef(*enu, 0, 90, -100) == approx(xyz)


def test_array_enu():
    np = pytest.importorskip("numpy")

    xyz = (np.asarray(0), np.asarray(A), np.asarray(50))
    llh = (np.asarray(0), np.asarray(90), np.asarray(-100))
    enu = pm.ecef2enu(*xyz, *llh)
    assert pm.enu2ecef(*enu, *llh) == approx(xyz)

    xyz = (np.atleast_1d(0), np.atleast_1d(A), np.atleast_1d(50))
    llh = (np.atleast_1d(0), np.atleast_1d(90), np.atleast_1d(-100))
    enu = pm.ecef2enu(*xyz, *llh)
    assert pm.enu2ecef(*enu, *llh) == approx(xyz)


@pytest.mark.parametrize(
    "enu,lla,xyz", [((0, 0, 0), (0, 0, 0), (A, 0, 0)), ((0, 0, 1000), (0, 0, 0), (A + 1000, 0, 0))]
)
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
