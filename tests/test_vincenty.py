#!/usr/bin/env python
import pytest
from pytest import approx
from pymap3d.vincenty import vreckon, vdist, track2
import numpy as np

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


def test_track2():
    lats, lons = track2(40, 80, 65, -148, npts=3)
    assert lats == approx([40, 69.633139886, 65])
    assert lons == approx([80, 113.06849104, -148])


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

    assert np.all(sr[0] == approx(asr))  # for older pytest
    assert aaz == approx(az)


if __name__ == '__main__':
    pytest.main(['-xrsv', __file__])
