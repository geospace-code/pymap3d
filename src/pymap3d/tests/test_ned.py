import pymap3d as pm
from pytest import approx

lla0 = (42, -82, 200)
aer0 = (33, 70, 1000)

ELL = pm.Ellipsoid.from_name("wgs84")
A = ELL.semimajor_axis
B = ELL.semiminor_axis


def test_ecef_ned():
    enu = pm.aer2enu(*aer0)
    ned = (enu[1], enu[0], -enu[2])
    xyz = pm.aer2ecef(*aer0, *lla0)

    ned1 = pm.ecef2ned(*xyz, *lla0)
    assert ned1 == approx(ned)

    assert pm.ned2ecef(*ned, *lla0) == approx(xyz)


def test_enuv_nedv():
    vx, vy, vz = (5, 3, 2)
    ve, vn, vu = (5.368859646588048, 3.008520763668120, -0.352347711524077)
    assert pm.ecef2enuv(vx, vy, vz, *lla0[:2]) == approx((ve, vn, vu))

    assert pm.ecef2nedv(vx, vy, vz, *lla0[:2]) == approx((vn, ve, -vu))


def test_ned_geodetic():
    lla1 = pm.aer2geodetic(*aer0, *lla0)

    enu3 = pm.geodetic2enu(*lla1, *lla0)
    ned3 = (enu3[1], enu3[0], -enu3[2])

    assert pm.geodetic2ned(*lla1, *lla0) == approx(ned3)

    lla2 = pm.enu2geodetic(*enu3, *lla0)
    assert lla2 == approx(lla1)
    assert all(isinstance(n, float) for n in lla2)

    lla2 = pm.ned2geodetic(*ned3, *lla0)
    assert lla2 == approx(lla1)
    assert all(isinstance(n, float) for n in lla2)
