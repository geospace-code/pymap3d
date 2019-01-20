#!/usr/bin/env python
# Copyright (c) 2014-2018 Michael Hirsch, Ph.D.

"""
runs tests
"""
import pytest
from pytest import approx
import numpy as np
import pymap3d as pm

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

E = pm.Ellipsoid()


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


def test_aer_ecef():
    xyz = pm.aer2ecef(*aer0, *lla0)

    assert xyz == approx(axyz0)
    assert pm.aer2ecef(*raer0, *rlla0, deg=False) == approx(axyz0)

    with pytest.raises(ValueError):
        pm.aer2ecef(aer0[0], aer0[1], -1, *lla0)

    assert pm.ecef2aer(*xyz, *lla0) == approx(aer0)
    assert pm.ecef2aer(*xyz, *rlla0, deg=False) == approx(raer0)


@pytest.mark.parametrize('xyz, lla, aer', [((E.a - 1, 0, 0), (0, 0, 0), (0, -90, 1)),
                                           ((-E.a + 1, 0, 0), (0, 180, 0), (0, -90, 1)),
                                           ((0, E.a - 1, 0), (0, 90, 0), (0, -90, 1)),
                                           ((0, -E.a + 1, 0), (0, -90, 0), (0, -90, 1)),
                                           ((0, 0, E.b - 1), (90, 0, 0), (0, -90, 1)),
                                           ((0, 0, -E.b + 1), (-90, 0, 0), (0, -90, 1)), ])
def test_ecef2aer(xyz, lla, aer):
    assert pm.ecef2aer(*xyz, *lla) == approx(aer)


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
