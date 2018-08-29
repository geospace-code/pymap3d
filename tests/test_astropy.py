#!/usr/bin/env python
import pytest
from pytest import approx
import numpy as np
import pymap3d as pm
import pymap3d.datetime2hourangle as pmd
import pymap3d.haversine as pmh

lon = -148
t0 = '2014-04-06T08:00:00'
lla0 = (42, -82, 200)
sra = 2.90658
ha = 45.482789587392013
eci0 = (-3.977913815668146e6, -2.582332196263046e6, 4.250818828152067e6)


def test_datetime2sidereal():
    pytest.importorskip('astropy')
    # http://www.jgiesen.de/astro/astroJS/siderealClock/
    assert pmd.datetime2sidereal(t0, np.radians(lon), False) == approx(sra, rel=1e-5)

    assert pmd.datetime2sidereal([t0], np.radians(lon), False) == approx([sra], rel=1e-5)


def test_datetime2sidereal_vallado():
    assert pmd.datetime2sidereal(t0, np.radians(lon), True) == approx(sra, rel=1e-5)

    assert pmd.datetime2sidereal([t0], np.radians(lon), True) == approx([sra], rel=1e-5)


def test_anglesep():
    pytest.importorskip('astropy')

    assert pmh.anglesep(35, 23, 84, 20) == approx(ha)


def test_anglesep_meeus():
    # %% compare with astropy
    assert pmh.anglesep_meeus(35, 23, 84, 20) == approx(ha)


def test_eci():
    pytest.importorskip('astropy')

    t = '2013-01-15T12:00:05'
    lla = pm.eci2geodetic(eci0, t)
    assert lla == approx(lla0, rel=0.2)

    eci1 = pm.eci2ecef(eci0, t)
    assert eci1 == approx([649012.04640917, -4697980.55129606, 4250818.82815207])

    assert pm.ecef2eci(eci1, t) == approx(eci0)

    aer1 = pm.eci2aer(eci0, 42, -100, 0, t)
    assert aer1 == approx([83.73050, -6.614478, 1.473510e6])

    assert pm.aer2eci(*aer1, 42, -100, 0, t) == approx(eci0)

    with pytest.raises(ValueError):
        pm.aer2eci(aer1[0], aer1[1], -1, 42, -100, 0, t)


def test_eci_times():
    pytest.importorskip('astropy')

    with pytest.raises(AssertionError):
        pm.eci2ecef(eci0, [t0, t0])

    with pytest.raises(AssertionError):
        pm.ecef2eci(eci0, [t0, t0])

    eci0s = np.stack((eci0, eci0))
    assert pm.ecef2eci(pm.eci2ecef(eci0s, [t0] * 2), [t0] * 2) == approx(eci0s)


if __name__ == '__main__':
    pytest.main(['-xrsv', __file__])
