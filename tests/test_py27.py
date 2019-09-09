"""
for strictly pure Python (no Numpy) only modules
Python >= 2.7
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
    xyz = pm.aer2ecef(aer0[0], aer0[1], aer0[2], lla0[0], lla0[1], lla0[2])

    enu = pm.aer2enu(aer0[0], aer0[1], aer0[2])

    assert pm.enu2ecef(enu[0], enu[1], enu[2], lla0[0], lla0[1], lla0[2]) == approx(xyz)


if __name__ == '__main__':
    pytest.main(["-v", __file__])
