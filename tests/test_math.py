"""
for strictly pure Python (no Numpy) only modules
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


if __name__ == '__main__':
    pytest.main([__file__])
