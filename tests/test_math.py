"""
for strictly pure Python (no Numpy) only modules
Python >= 3.5
"""
import pytest
from pytest import approx
from math import radians

import pymap3d.math as pm

lla0 = (42, -82, 200)
rlla0 = (radians(lla0[0]), radians(lla0[1]), lla0[2])
aer0 = (33, 70, 1000)
raer0 = (radians(aer0[0]), radians(aer0[1]), aer0[2])
enu0 = (186.277521, 286.842228, 939.692621)
ned0 = (enu0[1], enu0[0], -enu0[2])
vx, vy, vz = (5, 3, 2)
ve, vn, vu = (5.368859646588048, 3.008520763668120, -0.352347711524077)


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
    pytest.main([__file__])
