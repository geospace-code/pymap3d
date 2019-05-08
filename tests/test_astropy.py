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


@pytest.mark.parametrize('useastropy', [True, False])
def test_eci_geodetic(useastropy):
    t = '2013-01-15T12:00:05'
    lla = pm.eci2geodetic(*eci0, t, useastropy=useastropy)
    assert lla == approx(lla0, rel=0.2)


@pytest.mark.parametrize('useastropy', [True, False])
def test_eci_aer(useastropy):
    t = '2013-01-15T12:00:05'

    aer1 = pm.eci2aer(*eci0, 42, -100, 0, t, useastropy=useastropy)
    assert aer1 == approx([83.73050, -6.614478, 1.473510e6], rel=0.001)

    assert pm.aer2eci(*aer1, 42, -100, 0, t, useastropy=useastropy) == approx(eci0, rel=0.001)

    with pytest.raises(ValueError):
        pm.aer2eci(aer1[0], aer1[1], -1, 42, -100, 0, t, useastropy=useastropy)


if __name__ == '__main__':
    pytest.main(['-xrsv', __file__])
