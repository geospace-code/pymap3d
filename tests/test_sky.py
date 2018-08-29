import pytest
from pytest import approx
import pymap3d as pm

lat, lon = (65, -148)
lla0 = (42, -82, 200)
azi, eli = (180.1, 80)
t0 = '2014-04-06T08:00:00'
ra, dec = (166.5032081149338, 55.000011165405752)


def test_azel2radec():
    R, D = pm.azel2radec(azi, eli, lat, lon, t0)
    assert R == approx(ra, rel=0.01)
    assert D == approx(dec, rel=0.01)

    Rv, Dv = pm.azel2radec(azi, eli, lat, lon, t0, usevallado=True)
    assert Rv == approx(ra)
    assert Dv == approx(dec)


def test_radec2azel():
    azapy, elapy = pm.radec2azel(ra, dec, lat, lon, t0)
    assert azapy == approx(azi, rel=0.01, abs=0.1)
    assert elapy == approx(eli, rel=0.01, abs=0.1)

    azvallado, elvallado = pm.radec2azel(ra, dec, lat, lon, t0, usevallado=True)
    assert azvallado == approx(azi, rel=1e-2)
    assert elvallado == approx(eli, rel=1e-2)


if __name__ == '__main__':
    pytest.main(['-xrsv', __file__])
