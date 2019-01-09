#!/usr/bin/env python
import pytest
from pytest import approx
import pymap3d as pm
from math import radians, nan
import numpy as np

lla0 = (42, -82, 200)
rlla0 = (radians(lla0[0]), radians(lla0[1]), lla0[2])
lla1 = (42.002582, -81.997752, 1.1397018e3)
rlla1 = (np.radians(lla1[0]), np.radians(lla1[1]), lla1[2])

xyz0 = (660675.2518247,
        -4700948.68316,
        4245737.66222)

aer0 = (33, 70, 1000)
raer0 = (np.radians(aer0[0]), np.radians(aer0[1]), aer0[2])

E = pm.Ellipsoid()


def test_ecef():
    xyz = pm.geodetic2ecef(*lla0)

    assert xyz == approx(xyz0)
    assert pm.geodetic2ecef(*rlla0, deg=False) == approx(xyz)

    with pytest.raises(ValueError):
        pm.geodetic2ecef(-100, lla0[1], lla0[2])

    with pytest.raises(ValueError):
        pm.geodetic2ecef(lla0[0], -200, lla0[2])

    assert pm.ecef2geodetic(*xyz) == approx(lla0)
    assert pm.ecef2geodetic(*xyz, deg=False) == approx(rlla0)

    assert pm.ecef2geodetic((E.a - 1) / np.sqrt(2),
                            (E.a - 1) / np.sqrt(2), 0)


def test_aer():
    lla2 = pm.aer2geodetic(*aer0, *lla0)
    rlla2 = pm.aer2geodetic(*raer0, *rlla0, deg=False)

    with pytest.raises(ValueError):
        pm.aer2geodetic(aer0[0], aer0[1], -1, *lla0)

    assert lla2 == approx(lla1)
    assert rlla2 == approx(rlla1)

    assert pm.geodetic2aer(*lla2, *lla0) == approx(aer0)
    assert pm.geodetic2aer(*rlla2, *rlla0, deg=False) == approx(raer0)


def test_allnan():

    anan = np.empty((10, 10))
    anan.fill(nan)
    assert np.isnan(pm.geodetic2aer(anan, anan, anan, *lla0)).all()
    assert np.isnan(pm.aer2geodetic(anan, anan, anan, *lla0)).all()


def test_somenan():
    xyz = np.stack((xyz0, (nan, nan, nan)))

    lat, lon, alt = pm.ecef2geodetic(xyz[:, 0], xyz[:, 1], xyz[:, 2])
    assert (lat[0], lon[0], alt[0]) == approx(lla0)


if __name__ == '__main__':
    pytest.main([__file__])
