import pytest
from pytest import approx

try:
    from numpy import asarray
except ImportError:

    def asarray(*args):  # type: ignore
        "dummy function to convert values to arrays"
        return args


import pymap3d as pm

ELL = pm.Ellipsoid.from_name("wgs84")
A = ELL.semimajor_axis
B = ELL.semiminor_axis

llrlla = [
    ((0, 0, A - 1), (0, 0, -1)),
    ((0, 90, A - 1), (0, 90, -1)),
    ((0, -90, A + 1), (0, -90, 1)),
    ((44.807576814237606, 270, 6367490.543857), (45, 270, 1)),
    ((90, 0, B + 1), (90, 0, 1)),
    ((90, 15, B - 1), (90, 15, -1)),
    ((-90, 0, B + 1), (-90, 0, 1)),
]
llallr = [
    ((0, 0, -1), (0, 0, A - 1)),
    ((0, 90, -1), (0, 90, A - 1)),
    ((0, -90, 1), (0, -90, A + 1)),
    ((45, 270, 1), (44.807576814237606, 270, 6367490.543857)),
    ((90, 0, 1), (90, 0, B + 1)),
    ((90, 15, -1), (90, 15, B - 1)),
    ((-90, 0, 1), (-90, 0, B + 1)),
]
llallr_list = [([[i] for i in lla], llr) for lla, llr in llallr]
llrlla_list = [([[i] for i in llr], lla) for llr, lla in llrlla]
llallr_array = [([asarray(i) for i in lla], llr) for lla, llr in llallr]
llrlla_array = [([asarray(i) for i in llr], lla) for llr, lla in llrlla]

atol_dist = 1e-6  # 1 micrometer


@pytest.mark.parametrize("lla, llr", llallr)
def test_geodetic2spherical(lla, llr):
    coords = pm.geodetic2spherical(*lla)
    assert coords[:2] == approx(llr[:2])
    assert coords[2] == approx(llr[2], abs=atol_dist)


@pytest.mark.parametrize("llr, lla", llrlla)
def test_spherical2geodetic(llr, lla):
    coords = pm.spherical2geodetic(*llr)
    assert coords[:2] == approx(lla[:2])
    assert coords[2] == approx(lla[2], abs=atol_dist)


@pytest.mark.parametrize("lla, llr", llallr_list)
def test_geodetic2spherical_list(lla, llr):
    pytest.importorskip("numpy")
    coords = pm.geodetic2spherical(*lla)
    assert coords[:2] == approx(llr[:2])
    assert coords[2] == approx(llr[2], abs=atol_dist)


@pytest.mark.parametrize("llr, lla", llrlla_list)
def test_spherical2geodetic_list(llr, lla):
    pytest.importorskip("numpy")
    coords = pm.spherical2geodetic(*llr)
    assert coords[:2] == approx(lla[:2])
    assert coords[2] == approx(lla[2], abs=atol_dist)


@pytest.mark.parametrize("lla, llr", llallr_array)
def test_geodetic2spherical_array(lla, llr):
    pytest.importorskip("numpy")
    coords = pm.geodetic2spherical(*lla)
    assert coords[:2] == approx(llr[:2])
    assert coords[2] == approx(llr[2], abs=atol_dist)


@pytest.mark.parametrize("llr, lla", llrlla_array)
def test_spherical2geodetic_array(llr, lla):
    pytest.importorskip("numpy")
    coords = pm.spherical2geodetic(*llr)
    assert coords[:2] == approx(lla[:2])
    assert coords[2] == approx(lla[2], abs=atol_dist)
