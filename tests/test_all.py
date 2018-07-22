#!/usr/bin/env python
# Copyright (c) 2014-2018 Michael Hirsch, Ph.D.

"""
runs tests
"""
import pytest
from datetime import datetime
import numpy as np
from numpy.testing import assert_allclose
from pymap3d.vincenty import vreckon, vdist

try:  # for validation
    import pyproj
except ImportError:
    pyproj = None

import pymap3d as pm
from pymap3d.timeconv import str2dt

lla0 = (42, -82, 200)
rlla0 = (np.radians(lla0[0]), np.radians(lla0[1]), lla0[2])

aer0 = (33, 70, 1000)
raer0 = (np.radians(aer0[0]), np.radians(aer0[1]), aer0[2])

# %% outcomes from matlab
xyz0 = (660.6753e3, -4700.9487e3, 4245.738e3)  # geodetic2ecef

lla1 = (42.002582, -81.997752, 1.1397018e3)  # aer2geodetic
rlla1 = (np.radians(lla1[0]), np.radians(lla1[1]), lla1[2])

axyz0 = 660930.2, -4701424, 4246579.6  # aer2ecef

enu0 = (186.277521, 286.842228, 939.692621)  # aer2enu
ned0 = (enu0[1], enu0[0], -enu0[2])

# vector
vx, vy, vz = (5, 3, 2)
ve, vn, vu = (5.368859646588048, 3.008520763668120, -0.352347711524077)


def test_str2dt():

    assert str2dt(datetime(2014, 4, 6, 8)) == datetime(2014, 4, 6, 8)  # passthrough
    assert str2dt('2014-04-06T08:00:00') == datetime(2014, 4, 6, 8)
    ti = [str2dt('2014-04-06T08:00:00'), str2dt('2014-04-06T08:01:02')]
    to = [datetime(2014, 4, 6, 8), datetime(2014, 4, 6, 8, 1, 2)]
    assert ti == to   # even though ti is numpy array of datetime and to is list of datetime

# %%


def test_geodetic():
    if pyproj:
        ecef = pyproj.Proj(proj='geocent', ellps='WGS84', datum='WGS84')
        lla = pyproj.Proj(proj='latlong', ellps='WGS84', datum='WGS84')

    xyz1 = pm.geodetic2ecef(*lla0)

    assert_allclose(pm.geodetic2ecef(*rlla0, deg=False), xyz1, err_msg='geodetic2ecef: rad')
    assert_allclose(xyz1, xyz0, err_msg='geodetic2ecef: deg')

    assert_allclose(pm.ecef2geodetic(*xyz1), lla0, err_msg='ecef2geodetic: deg')
    assert_allclose(pm.ecef2geodetic(*xyz1, deg=False), rlla0, err_msg='ecef2geodetic: rad')

    if pyproj:
        assert_allclose(pyproj.transform(lla, ecef, lla0[1], lla0[0], lla0[2]), xyz1)
        assert_allclose(pyproj.transform(ecef, lla, *xyz1),
                        (lla0[1], lla0[0], lla0[2]))

    lla2 = pm.aer2geodetic(*aer0, *lla0)
    rlla2 = pm.aer2geodetic(*raer0, *rlla0, deg=False)

    assert_allclose(lla2, lla1, err_msg='aer2geodetic: deg')
    assert_allclose(rlla2, rlla1, err_msg='aer2geodetic:rad')

    assert_allclose(pm.geodetic2aer(*lla2, *lla0), aer0, err_msg='geodetic2aer: deg')
    assert_allclose(pm.geodetic2aer(*rlla2, *rlla0, deg=False), raer0, err_msg='geodetic2aer: rad')

# %% aer-ecef
    xyz2 = pm.aer2ecef(*aer0, *lla0)

    assert_allclose(pm.aer2ecef(*raer0, *rlla0, deg=False),
                    axyz0, err_msg='aer2ecef:rad')

    assert_allclose(xyz2, axyz0, err_msg='aer2ecef: deg')

    assert_allclose(pm.ecef2aer(*xyz2, *lla0), aer0, err_msg='ecef2aer:deg')
    assert_allclose(pm.ecef2aer(*xyz2, *rlla0, deg=False), raer0, err_msg='ecef2aer:rad')
# %% aer-enu
    enu1 = pm.aer2enu(*aer0)
    ned1 = (enu1[1], enu1[0], -enu1[2])

    assert_allclose(enu1, enu0, err_msg='aer2enu: deg')
    assert_allclose(pm.aer2enu(*raer0, deg=False), enu0, err_msg='aer2enu: rad')

    assert_allclose(pm.aer2ned(*aer0), ned0, err_msg='aer2ned')

    assert_allclose(pm.enu2aer(*enu1), aer0, err_msg='enu2aer: deg')
    assert_allclose(pm.enu2aer(*enu1, deg=False), raer0, err_msg='enu2aer: rad')

    assert_allclose(pm.ned2aer(*ned1), aer0, err_msg='ned2aer')

# %% enu-ecef
    assert_allclose(pm.enu2ecef(*enu1, *lla0), xyz2, err_msg='enu2ecef: deg')
    assert_allclose(pm.enu2ecef(*enu1, *rlla0, deg=False), xyz2, err_msg='enu2ecef: rad')

    assert_allclose(pm.ecef2enu(*xyz2, *lla0), enu1, err_msg='ecef2enu:deg')
    assert_allclose(pm.ecef2enu(*xyz2, *rlla0, deg=False), enu1, err_msg='ecef2enu:rad')

    assert_allclose(pm.ecef2ned(*xyz2, *lla0), ned1,
                    err_msg='ecef2ned')

    assert_allclose(pm.ned2ecef(*ned1, *lla0), xyz2,
                    err_msg='ned2ecef')
# %%
    assert_allclose(pm.ecef2enuv(vx, vy, vz, *lla0[:2]), (ve, vn, vu))

    assert_allclose(pm.ecef2nedv(vx, vy, vz, *lla0[:2]), (vn, ve, -vu))

# %%
    enu3 = pm.geodetic2enu(*lla2, *lla0)
    ned3 = (enu3[1], enu3[0], -enu3[2])

    assert_allclose(pm.geodetic2ned(*lla2, *lla0), ned3,
                    err_msg='geodetic2ned: deg')

    assert_allclose(pm.enu2geodetic(*enu3, *lla0), lla2,
                    err_msg='enu2geodetic')

    assert_allclose(pm.ned2geodetic(*ned3, *lla0), lla2,
                    err_msg='ned2geodetic')


def test_vincenty():
    az = 38
    sr = 3e3
    lat2, lon2, a21 = vreckon(10, 20, sr, az)
    assert_allclose((lat2, lon2, a21),
                    (10.02137267, 20.016847, 218.0029286))
    if pyproj:
        p4lon, p4lat, p4a21 = pyproj.Geod(ellps='WGS84').fwd(lon2, lat2, az, sr)
        assert_allclose((p4lon, p4lat, p4a21 % 360.), (lon2, lat2, a21), rtol=0.0025)

    assert_allclose(vdist(10, 20, lat2, lon2), (sr, az, a21))
    if pyproj:
        p4az, p4a21, p4sr = pyproj.Geod(ellps='WGS84').inv(20, 10, lon2, lat2)
        assert_allclose((p4az, p4a21 % 360., p4sr), (az, a21, sr))


if __name__ == '__main__':
    pytest.main()
