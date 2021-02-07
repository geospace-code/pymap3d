import pytest
from pytest import approx
import pymap3d as pm

ell = pm.Ellipsoid()
A = ell.semimajor_axis


@pytest.mark.parametrize(
    "lat,curvature", [(0, A), (90, 0), (-90, 0), (45.0, 4517590.87884893), (-45, 4517590.87884893)]
)
def test_rcurve_parallel(lat, curvature):
    assert pm.rcurve_parallel(lat) == approx(curvature, abs=1e-9)


def test_numpy_parallel():
    pytest.importorskip("numpy")
    assert pm.rcurve_parallel([0, 90]) == approx([A, 0], abs=1e-9)


@pytest.mark.parametrize(
    "lat,curvature",
    [
        (0, 6335439.327),
        (90, 6399593.6258),
        (-90, 6399593.6258),
        (45.0, 6367381.8156),
        (-45, 6367381.8156),
    ],
)
def test_rcurve_meridian(lat, curvature):
    assert pm.rcurve_meridian(lat) == approx(curvature)


def test_numpy_meridian():
    pytest.importorskip("numpy")
    assert pm.rcurve_meridian([0, 90]) == approx([6335439.327, 6399593.6258])


def test_numpy_transverse():
    pytest.importorskip("numpy")
    assert pm.rcurve_transverse([0, 90]) == approx([A, 6399593.626])
