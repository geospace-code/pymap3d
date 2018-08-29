#!/usr/bin/env python
# Copyright (c) 2014-2018 Michael Hirsch, Ph.D.

"""
runs tests
"""
import pytest
from pytest import approx
from datetime import datetime
import numpy as np
import pymap3d as pm
from pymap3d.timeconv import str2dt
from packaging import version

pi = np.pi
nan = np.nan

OLD = version.parse(pytest.__version__) < version.parse('3.5')

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


def test_ellipsoid():

    assert pm.ecef2geodetic(*xyz0, ell=pm.Ellipsoid('wgs84')) == approx([42., -82., 200.24339])
    assert pm.ecef2geodetic(*xyz0, ell=pm.Ellipsoid('grs80')) == approx([42., -82., 200.24344])
    assert pm.ecef2geodetic(*xyz0, ell=pm.Ellipsoid('clrk66')) == approx([42.00213, -82., 237.17182])
    assert pm.ecef2geodetic(*xyz0, ell=pm.Ellipsoid('mars')) == approx([41.99476, -82., 2.981169e6])
    assert pm.ecef2geodetic(*xyz0, ell=pm.Ellipsoid('venus')) == approx([41.808706, -82., 3.178069e5])
    assert pm.ecef2geodetic(*xyz0, ell=pm.Ellipsoid('moon')) == approx([41.808706, -82., 4.630807e6])


def test_str2dt():

    assert str2dt(datetime(2014, 4, 6, 8)) == datetime(2014, 4, 6, 8)  # passthrough
    assert str2dt('2014-04-06T08:00:00') == datetime(2014, 4, 6, 8)
    ti = [str2dt('2014-04-06T08:00:00'), str2dt('2014-04-06T08:01:02')]
    to = [datetime(2014, 4, 6, 8), datetime(2014, 4, 6, 8, 1, 2)]
    assert ti == to   # even though ti is numpy array of datetime and to is list of datetime

# %%


@pytest.mark.skipif(OLD, reason='too old pytest')
def test_losint():

    az = [0., 10., 125.]

    lat, lon, sr = pm.lookAtSpheroid(*lla0, az, tilt=0.)
    assert (lat[0] == lat).all() and (lon[0] == lon).all() and (sr[0] == sr).all()

    assert (lat[0], lon[0], sr[0]) == approx(lla0)

    with pytest.raises(ValueError):
        pm.lookAtSpheroid(lla0[0], lla0[1], -1, az, 0)

# %%
    tilt = [30., 45., 90.]
    lat, lon, sr = pm.lookAtSpheroid(*lla0, az, tilt)

    truth = np.array([[42.00103959, lla0[1], 230.9413173],
                      [42.00177328, -81.9995808, 282.84715651],
                      [nan, nan, nan]])

    assert np.column_stack((lat, lon, sr)) == approx(truth, nan_ok=True)


def test_geodetic():
    xyz = pm.geodetic2ecef(*lla0)

    assert xyz == approx(xyz0)
    assert pm.geodetic2ecef(*rlla0, deg=False) == approx(xyz)

    with pytest.raises(ValueError):
        pm.geodetic2ecef(lla0[0], lla0[1], -1)

    with pytest.raises(ValueError):
        pm.geodetic2ecef(-100, lla0[1], lla0[2])

    with pytest.raises(ValueError):
        pm.geodetic2ecef(lla0[0], -200, lla0[2])

    assert pm.ecef2geodetic(*xyz) == approx(lla0)
    assert pm.ecef2geodetic(*xyz, deg=False) == approx(rlla0)

    lla2 = pm.aer2geodetic(*aer0, *lla0)
    rlla2 = pm.aer2geodetic(*raer0, *rlla0, deg=False)

    with pytest.raises(ValueError):
        pm.aer2geodetic(aer0[0], aer0[1], -1, *lla0)

    assert lla2 == approx(lla1)
    assert rlla2 == approx(rlla1)

    assert pm.geodetic2aer(*lla2, *lla0) == approx(aer0)
    assert pm.geodetic2aer(*rlla2, *rlla0, deg=False) == approx(raer0)

    anan = np.empty((10, 10))
    anan.fill(np.nan)
    assert np.isnan(pm.geodetic2aer(anan, anan, anan, *lla0)).all()
    assert np.isnan(pm.aer2geodetic(anan, anan, anan, *lla0)).all()


def test_aer_ecef():
    xyz = pm.aer2ecef(*aer0, *lla0)

    assert xyz == approx(axyz0)
    assert pm.aer2ecef(*raer0, *rlla0, deg=False) == approx(axyz0)

    with pytest.raises(ValueError):
        pm.aer2ecef(aer0[0], aer0[1], -1, *lla0)

    assert pm.ecef2aer(*xyz, *lla0) == approx(aer0)
    assert pm.ecef2aer(*xyz, *rlla0, deg=False) == approx(raer0)


def test_aer_enu():
    xyz = pm.aer2ecef(*aer0, *lla0)

    enu = pm.aer2enu(*aer0)

    assert enu == approx(enu0)
    assert pm.aer2enu(*raer0, deg=False) == approx(enu0)

    with pytest.raises(ValueError):
        pm.aer2enu(aer0[0], aer0[1], -1)

    assert pm.enu2ecef(*enu, *lla0) == approx(xyz)
    assert pm.enu2ecef(*enu, *rlla0, deg=False) == approx(xyz)

    assert pm.ecef2enu(*xyz, *lla0) == approx(enu)
    assert pm.ecef2enu(*xyz, *rlla0, deg=False) == approx(enu)


def test_ned():
    xyz = pm.aer2ecef(*aer0, *lla0)
    enu = pm.aer2enu(*aer0)
    ned = (enu[1], enu[0], -enu[2])
    lla = pm.aer2geodetic(*aer0, *lla0)

    assert pm.aer2ned(*aer0) == approx(ned0)

    with pytest.raises(ValueError):
        pm.aer2ned(aer0[0], aer0[1], -1)

    assert pm.enu2aer(*enu) == approx(aer0)
    assert pm.enu2aer(*enu, deg=False) == approx(raer0)

    assert pm.ned2aer(*ned) == approx(aer0)

    assert pm.ecef2ned(*xyz, *lla0) == approx(ned)

    assert pm.ned2ecef(*ned, *lla0) == approx(xyz)
# %%
    assert pm.ecef2enuv(vx, vy, vz, *lla0[:2]) == approx((ve, vn, vu))

    assert pm.ecef2nedv(vx, vy, vz, *lla0[:2]) == approx((vn, ve, -vu))

# %%
    enu3 = pm.geodetic2enu(*lla, *lla0)
    ned3 = (enu3[1], enu3[0], -enu3[2])

    assert pm.geodetic2ned(*lla, *lla0) == approx(ned3)

    assert pm.enu2geodetic(*enu3, *lla0) == approx(lla)

    assert pm.ned2geodetic(*ned3, *lla0) == approx(lla)


if __name__ == '__main__':
    pytest.main(['-xrsv', __file__])
