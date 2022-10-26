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

    n, e, d = pm.ecef2ned(*xyz, *lla0)
    assert n == approx(ned[0])
    assert e == approx(ned[1])
    assert d == approx(ned[2])

    assert pm.ned2ecef(*ned, *lla0) == approx(xyz)


def test_enuv_nedv():
    vx, vy, vz = (5, 3, 2)
    ve, vn, vu = (5.368859646588048, 3.008520763668120, -0.352347711524077)
    assert pm.ecef2enuv(vx, vy, vz, *lla0[:2]) == approx((ve, vn, vu))

    assert pm.ecef2nedv(vx, vy, vz, *lla0[:2]) == approx((vn, ve, -vu))


def test_ned_geodetic():
    lat1, lon1, alt1 = pm.aer2geodetic(*aer0, *lla0)

    enu3 = pm.geodetic2enu(lat1, lon1, alt1, *lla0)
    ned3 = (enu3[1], enu3[0], -enu3[2])

    assert pm.geodetic2ned(lat1, lon1, alt1, *lla0) == approx(ned3)

    lat, lon, alt = pm.enu2geodetic(*enu3, *lla0)
    assert lat == approx(lat1)
    assert lon == approx(lon1)
    assert alt == approx(alt1)
    assert isinstance(lat, float)
    assert isinstance(lon, float)
    assert isinstance(alt, float)

    lat, lon, alt = pm.ned2geodetic(*ned3, *lla0)
    assert lat == approx(lat1)
    assert lon == approx(lon1)
    assert alt == approx(alt1)
    assert isinstance(lat, float)
    assert isinstance(lon, float)
    assert isinstance(alt, float)
