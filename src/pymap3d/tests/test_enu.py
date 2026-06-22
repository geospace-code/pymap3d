from math import radians

import pymap3d as pm
import pytest
from pytest import approx


def get_ellipsoid_params():
    ell = pm.Ellipsoid.from_name("wgs84")
    return ell.semimajor_axis, ell.semiminor_axis


A, B = get_ellipsoid_params()


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

    x = np.asarray(0)
    y = np.asarray(A)
    z = np.asarray(50)
    lat = np.asarray(0)
    lon = np.asarray(90)
    h = np.asarray(-100)

    e, n, u = pm.ecef2enu(x, y, z, lat, lon, h)

    assert pm.enu2ecef(e, n, u, lat, lon, h) == approx((x, y, z))

    x = np.atleast_1d(0)
    y = np.atleast_1d(A)
    z = np.atleast_1d(50)
    lat = np.atleast_1d(0)
    lon = np.atleast_1d(90)
    h = np.atleast_1d(-100)

    e, n, u = pm.ecef2enu(x, y, z, lat, lon, h)
    assert pm.enu2ecef(e, n, u, lat, lon, h) == approx((x, y, z))


@pytest.mark.parametrize(
    "enu,lla,xyz",
    [((0, 0, 0), (0, 0, 0), (A, 0, 0)), ((0, 0, 1000), (0, 0, 0), (A + 1000, 0, 0))],
)
def test_enu_ecef(enu, lla, xyz):
    x, y, z = pm.enu2ecef(*enu, *lla)
    assert x == approx(xyz[0])
    assert y == approx(xyz[1])
    assert z == approx(xyz[2])
    assert isinstance(x, float)
    assert isinstance(y, float)
    assert isinstance(z, float)

    rlat = radians(lla[0])
    rlon = radians(lla[1])
    rh = lla[2]

    assert pm.enu2ecef(*enu, rlat, rlon, rh, deg=False) == approx(xyz)

    e, n, u = pm.ecef2enu(*xyz, *lla)
    assert e == approx(enu[0])
    assert n == approx(enu[1])
    assert u == approx(enu[2])
    assert isinstance(e, float)
    assert isinstance(n, float)
    assert isinstance(u, float)

    e, n, u = pm.ecef2enu(*xyz, rlat, rlon, rh, deg=False)
    assert e == approx(enu[0])
    assert n == approx(enu[1])
    assert u == approx(enu[2])
