#!/usr/bin/env python
import pytest
from pytest import approx
from pymap3d.vincenty import vreckon, vdist
import pymap3d as pm

az = [38., 45.]
sr = [3e3, 1e3]

lla0 = [42, -82, 200]
ll0 = [10, 20]

lla1 = [40, -80, 1120]


lat2 = [10.02137267, 10.01917819]
lon2 = [20.0168471, 20.0193493]
az2 = [218.00292856, 225.00336316]

lat3 = (10.02137267, 10.00639286)
lon3 = (20.0168471, 20.00644951)
az3 = (218.00292856, 225.0011203)


def test_vreckon():
    """ tests scalars, vectors"""

    # scalar
    assert vreckon(*ll0, sr[0], az[0]) == approx((lat2[0], lon2[0], az2[0]))
    # az vector
    a, b, c = vreckon(*ll0, sr[0], az)
    assert a == approx(lat2)
    assert b == approx(lon2)
    assert c == approx(az2)
    # rng, az vectors
    a, b, c = vreckon(*ll0, sr, az)
    assert a == approx(lat3)
    assert b == approx(lon3)
    assert c == approx(az3)


def test_vdist():
    lat1, lon1, a21 = vreckon(*ll0, sr[0], az[0])

    assert vdist(*ll0, lat1, lon1) == approx((sr[0], az[0], a21))
    # lat, lon vectors
    asr, aaz, aa21 = vdist(*ll0, lat2, lon2)

    assert asr == approx(sr[0])
    assert aaz == approx(az)


def test_compare_vicenty():
    taz, tsr = az[0], sr[0]
    pyproj = pytest.importorskip('pyproj')

    lat2, lon2, a21 = vreckon(10, 20, tsr, taz)

    p4lon, p4lat, p4a21 = pyproj.Geod(ellps='WGS84').fwd(lon2, lat2, taz, tsr)
    assert (p4lon, p4lat, p4a21 % 360.) == approx((lon2, lat2, a21), rel=0.0025)

    p4az, p4a21, p4sr = pyproj.Geod(ellps='WGS84').inv(20, 10, lon2, lat2)
    assert (p4az, p4a21 % 360., p4sr) == approx((taz, a21, tsr))


def test_compare_geodetic():
    pyproj = pytest.importorskip('pyproj')

    xyz = pm.geodetic2ecef(*lla0)

    ecef = pyproj.Proj(proj='geocent', ellps='WGS84', datum='WGS84')
    lla = pyproj.Proj(proj='latlong', ellps='WGS84', datum='WGS84')

    assert pyproj.transform(lla, ecef, lla0[1], lla0[0], lla0[2]) == approx(xyz)
    assert pyproj.transform(ecef, lla, *xyz) == approx((lla0[1], lla0[0], lla0[2]))


if __name__ == '__main__':
    pytest.main(['-xv', __file__])
