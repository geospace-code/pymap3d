import pytest
from pytest import approx
import numpy as np
from datetime import datetime
import pymap3d as pm
import pymap3d.datetime2hourangle as pmd
import pymap3d.haversine as pmh
import pymap3d.vallado as pmv
from numpy.testing import assert_allclose

t0 = '2014-04-06T08:00:00'
lat, lon = (65, -148)
lla0 = (42, -82, 200)
azi, eli = (180.1, 80)
ra, dec = (166.5032081149338, 55.000011165405752)
sra = 2.90658
ha = 45.482789587392013


def test_datetime2sidereal():
    pytest.importorskip('astropy')
    # http://www.jgiesen.de/astro/astroJS/siderealClock/
    assert_allclose(pmd.datetime2sidereal(t0, np.radians(lon), False), sra, rtol=1e-5)

    assert_allclose(pmd.datetime2sidereal([t0], np.radians(lon), False), [sra], rtol=1e-5)


def test_datetime2sidereal_vallado():
    assert_allclose(pmd.datetime2sidereal(t0, np.radians(lon), True), sra, rtol=1e-5)

    assert_allclose(pmd.datetime2sidereal([t0], np.radians(lon), True), [sra], rtol=1e-5)


def test_anglesep():
    pytest.importorskip('astropy')

    assert pmh.anglesep(35, 23, 84, 20) == approx(ha)


def test_anglesep_meeus():
    # %% compare with astropy
    assert_allclose(pmh.anglesep_meeus(35, 23, 84, 20), ha)


def test_eci():
    pytest.importorskip('astropy')

    teci = (-3.977913815668146e6, -2.582332196263046e6, 4.250818828152067e6)
    t = datetime(2013, 1, 15, 12, 0, 5)
    lla = np.asarray(pm.eci2geodetic(teci, t)).squeeze()
    assert_allclose(lla, lla0, rtol=0.2)

    assert_allclose(pm.eci2ecef(teci, t).squeeze(),
                    [649012.04640917, -4697980.55129606, 4250818.82815207])

    assert pm.ecef2eci([649012.04640917, -4697980.55129606, 4250818.82815207], t).squeeze() == approx(teci)

    assert_allclose(np.asarray(pm.eci2aer(teci, 42, -100, 0, t)).squeeze(),
                    [83.73050, -6.614478, 1.473510e6])


def test_azel2radec():
    R, D = pm.azel2radec(azi, eli, lat, lon, t0)
    assert_allclose(R, ra, rtol=1e-2)
    assert_allclose(D, dec, rtol=1e-2)

    Rv, Dv = pmv.vazel2radec(azi, eli, lat, lon, t0)
    assert Rv == approx(ra)
    assert Dv == approx(dec)


def test_radec2azel():
    azapy, elapy = pm.radec2azel(ra, dec, lat, lon, t0)
    assert azapy == approx(azi, rel=0.01, abs=0.1)
    assert elapy == approx(eli, rel=0.01, abs=0.1)

    azvallado, elvallado = pmv.vradec2azel(ra, dec, lat, lon, t0)
    assert azvallado == approx(azi, rel=1e-2)
    assert elvallado == approx(eli, rel=1e-2)


if __name__ == '__main__':
    pytest.main(['-x', __file__])
