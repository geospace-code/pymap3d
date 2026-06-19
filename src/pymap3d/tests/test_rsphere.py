import pymap3d as pm
import pymap3d.rcurve as rcurve
import pymap3d.rsphere as rsphere
import pytest
from pytest import approx

ell = pm.Ellipsoid.from_name("wgs84")
A = ell.semimajor_axis


def test_geocentric_radius():
    assert rcurve.geocentric_radius(0) == approx(ell.semimajor_axis)
    assert rcurve.geocentric_radius(90) == approx(ell.semiminor_axis)
    assert rcurve.geocentric_radius(45) == approx(6367490.0)
    assert rcurve.geocentric_radius(30) == approx(6372824.0)


def test_rsphere_eqavol():
    assert rsphere.eqavol() == approx(6371000.8049)


def test_rsphere_authalic():
    assert rsphere.authalic() == approx(6371007.1809)

def test_rsphere_curve():
    assert rsphere.curve(0) == approx(6356788.163646411)
    assert rsphere.curve(90) == approx(6399593.625758492)
    assert rsphere.curve(45) == approx(6378110.052870348)
    assert rsphere.curve(30) == approx(6367429.010702812)


def test_rsphere_rectifying():
    assert rsphere.rectifying() == approx(6367449.1458)


def test_rsphere_biaxial():
    assert rsphere.biaxial() == approx(6367444.657)


def test_rsphere_triaxial():
    assert rsphere.triaxial() == approx(6371008.77)


def test_rsphere_euler():
    assert rsphere.euler(42, 82, 44, 100) == approx(6386606.829131)


def test_numpy_rsphere_euler():
    np = pytest.importorskip("numpy")
    lat1 = np.array([42, 0])
    lon1 = np.array([82, 0])
    assert rsphere.euler(lat1, lon1, 44, 100) == approx([6386606.829131, 6363111.70923164])
