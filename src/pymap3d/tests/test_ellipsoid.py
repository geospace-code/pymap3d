import pymap3d as pm
import pytest
from pytest import approx

xyz0 = (660e3, -4700e3, 4247e3)


@pytest.mark.parametrize(
    "model,f",
    [
        ("wgs84", 3.352810664747480e-03),
        ("wgs84_mean", 0.0),
        ("wgs72", 3.352779454167505e-03),
        ("grs80", 3.352810681182319e-03),
        ("clarke1866", 3.390075303928791e-03),
        ("moon", 0.0),
    ],
)
def test_reference(model, f):
    assert pm.Ellipsoid(model).flattening == approx(f)


def test_ellipsoid():

    assert pm.ecef2geodetic(*xyz0, ell=pm.Ellipsoid("wgs84")) == approx(
        [42.014670535, -82.0064785, 276.9136916]
    )
    assert pm.ecef2geodetic(*xyz0, ell=pm.Ellipsoid("wgs84_mean")) == approx(
        [41.823366301, -82.0064785, -2.13061272e3]
    )
    assert pm.ecef2geodetic(*xyz0, ell=pm.Ellipsoid("grs80")) == approx(
        [42.014670536, -82.0064785, 276.9137385]
    )
    assert pm.ecef2geodetic(*xyz0, ell=pm.Ellipsoid("clarke1866")) == approx(
        [42.01680003, -82.0064785, 313.9026793]
    )
    assert pm.ecef2geodetic(*xyz0, ell=pm.Ellipsoid("mars")) == approx(
        [42.009428417, -82.006479, 2.981246e6]
    )
    assert pm.ecef2geodetic(*xyz0, ell=pm.Ellipsoid("venus")) == approx(
        [41.8233663, -82.0064785, 3.17878159e5]
    )
    assert pm.ecef2geodetic(*xyz0, ell=pm.Ellipsoid("moon")) == approx(
        [41.8233663, -82.0064785, 4.630878e6]
    )
