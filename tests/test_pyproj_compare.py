#!/usr/bin/env python
import pytest
from pytest import approx
from pymap3d.vincenty import vreckon, vdist
import pymap3d as pm

az = 38.
sr = 3e3
lla0 = [42, -82, 200]
ll0 = [10, 20]

rlla2 = (10.02137267, 20.016847, 218.0029286)


def test_vincenty():
    """ tests scalars and vectors"""

    lat2, lon2, a21 = vreckon(*ll0, sr, az)
    assert (lat2, lon2, a21) == approx(rlla2)
    assert vreckon(*ll0, sr, [az, az]) == approx(rlla2)
    assert vreckon(*ll0, [sr, sr], [az, az]) == approx(rlla2)

    assert vdist(*ll0, lat2, lon2) == approx((sr, az, a21))
    assert vdist(*ll0, [lat2, lat2], [lon2, lon2]) == approx((sr, az, a21))


def test_compare_vicenty():
    pyproj = pytest.importorskip('pyproj')

    lat2, lon2, a21 = vreckon(10, 20, sr, az)

    p4lon, p4lat, p4a21 = pyproj.Geod(ellps='WGS84').fwd(lon2, lat2, az, sr)
    assert (p4lon, p4lat, p4a21 % 360.) == approx((lon2, lat2, a21), rel=0.0025)

    p4az, p4a21, p4sr = pyproj.Geod(ellps='WGS84').inv(20, 10, lon2, lat2)
    assert (p4az, p4a21 % 360., p4sr) == approx((az, a21, sr))


def test_compare_geodetic():
    pyproj = pytest.importorskip('pyproj')

    xyz = pm.geodetic2ecef(*lla0)

    ecef = pyproj.Proj(proj='geocent', ellps='WGS84', datum='WGS84')
    lla = pyproj.Proj(proj='latlong', ellps='WGS84', datum='WGS84')

    assert pyproj.transform(lla, ecef, lla0[1], lla0[0], lla0[2]) == approx(xyz)
    assert pyproj.transform(ecef, lla, *xyz) == approx((lla0[1], lla0[0], lla0[2]))


if __name__ == '__main__':
    pytest.main(['-xv', __file__])
