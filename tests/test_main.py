#!/usr/bin/env python
# Copyright (c) 2014-2018 Michael Hirsch, Ph.D.

"""
runs tests
"""
import pytest
from datetime import datetime
import numpy as np
from numpy.testing import assert_allclose
import pymap3d as pm
from pymap3d.timeconv import str2dt

pi = np.pi
nan = np.nan

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


def test_losint():

    az = [0., 10., 125.]

    lat, lon, sr = pm.lookAtSpheroid(*lla0, az, tilt=0.)
    assert (lat[0] == lat).all() and (lon[0] == lon).all() and (sr[0] == sr).all()

    assert_allclose((lat[0], lon[0], sr[0]), lla0, err_msg='los identity')
# %%
    tilt = [30., 45., 90.]
    lat, lon, sr = pm.lookAtSpheroid(*lla0, az, tilt)

    truth = np.array([[42.00103959, lla0[1], 230.9413173],
                      [42.00177328, -81.9995808, 282.84715651],
                      [nan, nan, nan]])

    assert_allclose(np.column_stack((lat, lon, sr)), truth)


def test_geodetic():
    xyz = pm.geodetic2ecef(*lla0)

    assert_allclose(pm.geodetic2ecef(*rlla0, deg=False), xyz, err_msg='geodetic2ecef: rad')
    assert_allclose(xyz, xyz0, err_msg='geodetic2ecef: deg')

    assert_allclose(pm.ecef2geodetic(*xyz), lla0, err_msg='ecef2geodetic: deg')
    assert_allclose(pm.ecef2geodetic(*xyz, deg=False), rlla0, err_msg='ecef2geodetic: rad')

    lla2 = pm.aer2geodetic(*aer0, *lla0)
    rlla2 = pm.aer2geodetic(*raer0, *rlla0, deg=False)

    assert_allclose(lla2, lla1, err_msg='aer2geodetic: deg')
    assert_allclose(rlla2, rlla1, err_msg='aer2geodetic:rad')

    assert_allclose(pm.geodetic2aer(*lla2, *lla0), aer0, err_msg='geodetic2aer: deg')
    assert_allclose(pm.geodetic2aer(*rlla2, *rlla0, deg=False), raer0, err_msg='geodetic2aer: rad')


def test_aer_ecef():
    xyz = pm.aer2ecef(*aer0, *lla0)

    assert_allclose(pm.aer2ecef(*raer0, *rlla0, deg=False),
                    axyz0, err_msg='aer2ecef:rad')

    assert_allclose(xyz, axyz0, err_msg='aer2ecef: deg')

    assert_allclose(pm.ecef2aer(*xyz, *lla0), aer0, err_msg='ecef2aer:deg')
    assert_allclose(pm.ecef2aer(*xyz, *rlla0, deg=False), raer0, err_msg='ecef2aer:rad')


def test_aer_enu():
    xyz = pm.aer2ecef(*aer0, *lla0)

    enu = pm.aer2enu(*aer0)

    assert_allclose(enu, enu0, err_msg='aer2enu: deg')
    assert_allclose(pm.aer2enu(*raer0, deg=False), enu0, err_msg='aer2enu: rad')

    assert_allclose(pm.enu2ecef(*enu, *lla0), xyz, err_msg='enu2ecef: deg')
    assert_allclose(pm.enu2ecef(*enu, *rlla0, deg=False), xyz, err_msg='enu2ecef: rad')

    assert_allclose(pm.ecef2enu(*xyz, *lla0), enu, err_msg='ecef2enu:deg')
    assert_allclose(pm.ecef2enu(*xyz, *rlla0, deg=False), enu, err_msg='ecef2enu:rad')


def test_ned():
    xyz = pm.aer2ecef(*aer0, *lla0)
    enu = pm.aer2enu(*aer0)
    ned = (enu[1], enu[0], -enu[2])
    lla = pm.aer2geodetic(*aer0, *lla0)

    assert_allclose(pm.aer2ned(*aer0), ned0, err_msg='aer2ned')

    assert_allclose(pm.enu2aer(*enu), aer0, err_msg='enu2aer: deg')
    assert_allclose(pm.enu2aer(*enu, deg=False), raer0, err_msg='enu2aer: rad')

    assert_allclose(pm.ned2aer(*ned), aer0, err_msg='ned2aer')

    assert_allclose(pm.ecef2ned(*xyz, *lla0), ned, err_msg='ecef2ned')

    assert_allclose(pm.ned2ecef(*ned, *lla0), xyz, err_msg='ned2ecef')
# %%
    assert_allclose(pm.ecef2enuv(vx, vy, vz, *lla0[:2]), (ve, vn, vu))

    assert_allclose(pm.ecef2nedv(vx, vy, vz, *lla0[:2]), (vn, ve, -vu))

# %%
    enu3 = pm.geodetic2enu(*lla, *lla0)
    ned3 = (enu3[1], enu3[0], -enu3[2])

    assert_allclose(pm.geodetic2ned(*lla, *lla0), ned3, err_msg='geodetic2ned: deg')

    assert_allclose(pm.enu2geodetic(*enu3, *lla0), lla, err_msg='enu2geodetic')

    assert_allclose(pm.ned2geodetic(*ned3, *lla0), lla, err_msg='ned2geodetic')


if __name__ == '__main__':
    pytest.main()
