#!/usr/bin/env python
import pytest
from pytest import approx
from math import radians
from datetime import datetime
import pymap3d as pm
import pymap3d.sidereal as pmd
import pymap3d.haversine as pmh

lon = -148
t0 = datetime(2014, 4, 6, 8)
lla0 = (42, -82, 200)
sra = 2.90658
ha = 45.482789587392013
eci0 = (-3.977913815668146e6, -2.582332196263046e6, 4.250818828152067e6)


@pytest.mark.parametrize("time", [t0, [t0]], ids=("scalar", "list"))
def test_sidereal(time):
    pytest.importorskip("astropy")
    # http://www.jgiesen.de/astro/astroJS/siderealClock/
    tsr = pmd.datetime2sidereal(time, radians(lon), False)
    if isinstance(tsr, list):
        tsr = tsr[0]
    assert tsr == approx(sra, rel=1e-5)


@pytest.mark.parametrize("time", [t0, [t0]], ids=("scalar", "list"))
def test_sidereal_vallado(time):
    tsr = pmd.datetime2sidereal(time, radians(lon), True)
    if isinstance(tsr, list):
        tsr = tsr[0]
    assert tsr == approx(sra, rel=1e-5)


def test_anglesep():
    pytest.importorskip("astropy")

    assert pmh.anglesep(35, 23, 84, 20) == approx(ha)


def test_anglesep_meeus():
    # %% compare with astropy
    assert pmh.anglesep_meeus(35, 23, 84, 20) == approx(ha)


def test_eci_geodetic():
    pytest.importorskip("astropy")
    t = "2013-01-15T12:00:05"
    lla = pm.eci2geodetic(*eci0, t)
    assert lla == approx(lla0, rel=0.001)


def test_eci_aer():
    pytest.importorskip("astropy")
    t = "2013-01-15T12:00:05"

    aer1 = pm.eci2aer(*eci0, 42, -100, 0, t)
    assert aer1 == approx([83.9511, -6.66787, 1485123.89], rel=0.001)

    assert pm.aer2eci(*aer1, 42, -100, 0, t) == approx(eci0, rel=0.001)

    with pytest.raises(ValueError):
        pm.aer2eci(aer1[0], aer1[1], -1, 42, -100, 0, t)


if __name__ == "__main__":
    pytest.main([__file__])
