#!/usr/bin/env python
import pytest
from pytest import approx
import numpy as np
import pymap3d as pm
import pymap3d.sidereal as pmd
import pymap3d.haversine as pmh

lon = -148
t0 = '2014-04-06T08:00:00'
lla0 = (42, -82, 200)
sra = 2.90658
ha = 45.482789587392013
eci0 = (-3.977913815668146e6, -2.582332196263046e6, 4.250818828152067e6)


@pytest.mark.parametrize('time', [t0, [t0]], ids=('scalar', 'list'))
def test_sidereal(time):
    pytest.importorskip('astropy')
    # http://www.jgiesen.de/astro/astroJS/siderealClock/
    assert pmd.datetime2sidereal(time, np.radians(lon), False) == approx(sra, rel=1e-5)


@pytest.mark.parametrize('time', [t0, [t0]], ids=('scalar', 'list'))
def test_sidereal_vallado(time):
    assert pmd.datetime2sidereal(time, np.radians(lon), True) == approx(sra, rel=1e-5)


def test_anglesep():
    pytest.importorskip('astropy')

    assert pmh.anglesep(35, 23, 84, 20) == approx(ha)


def test_anglesep_meeus():
    # %% compare with astropy
    assert pmh.anglesep_meeus(35, 23, 84, 20) == approx(ha)


def test_eci_astropy():
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


def test_eci_times_astropy():
    pytest.importorskip('astropy')

    with pytest.raises(AssertionError):
        pm.eci2ecef(eci0, [t0, t0])

    with pytest.raises(AssertionError):
        pm.ecef2eci(eci0, [t0, t0])

    eci0s = np.stack((eci0, eci0))
    assert pm.ecef2eci(pm.eci2ecef(eci0s, [t0] * 2), [t0] * 2) == approx(eci0s)


def test_eci_vallado():
    t = '2013-01-15T12:00:05'
    lla = pm.eci2geodetic(eci0, t, useastropy=False)
    assert lla == approx(lla0, rel=0.2)

    eci1 = pm.eci2ecef(eci0, t, useastropy=False)
    assert eci1 == approx([649012.04640917, -4697980.55129606, 4250818.82815207], rel=0.001)

    assert pm.ecef2eci(eci1, t, useastropy=False) == approx(eci0, rel=0.001)

    aer1 = pm.eci2aer(eci0, 42, -100, 0, t, useastropy=False)
    assert aer1 == approx([83.73050, -6.614478, 1.473510e6], rel=0.001)

    assert pm.aer2eci(*aer1, 42, -100, 0, t, useastropy=False) == approx(eci0, rel=0.001)

    with pytest.raises(ValueError):
        pm.aer2eci(aer1[0], aer1[1], -1, 42, -100, 0, t, useastropy=False)


def test_eci_times_vallado():
    with pytest.raises(AssertionError):
        pm.eci2ecef(eci0, [t0, t0], useastropy=False)

    with pytest.raises(AssertionError):
        pm.ecef2eci(eci0, [t0, t0], useastropy=False)

    eci0s = np.stack((eci0, eci0))
    assert pm.ecef2eci(pm.eci2ecef(eci0s, [t0] * 2, useastropy=False),
                       [t0] * 2, useastropy=False) == approx(eci0s, rel=0.001)


if __name__ == '__main__':
    pytest.main(['-xrsv', __file__])
