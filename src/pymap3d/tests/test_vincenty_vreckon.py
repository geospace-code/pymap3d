from math import radians

import pymap3d.vincenty as vincenty

import pytest
from pytest import approx

ll0 = [10, 20]
lat2 = [10.02137267, 10.01917819]
lon2 = [20.0168471, 20.0193493]
az2 = [218.00292856, 225.00336316]

sr1 = [3e3, 1e3]
az1 = [38, 45]
lat3 = (10.02137267, 10.00639286)
lon3 = (20.0168471, 20.00644951)
az3 = (218.00292856, 225.0011203)


@pytest.mark.parametrize("deg", [True, False])
@pytest.mark.parametrize(
    "lat,lon,srange,az,lato,lono",
    [
        (0, 0, 0, 0, 0, 0),
        (0, 0, 1.001875e7, 90, 0, 90),
        (0, 0, 1.001875e7, 270, 0, -90),
        (0, 0, 1.001875e7, -90, 0, -90),
        (0, 0, 2.00375e7, 90, 0, 180),
        (0, 0, 2.00375e7, 270, 0, -180),
        (0, 0, 2.00375e7, -90, 0, -180),
    ],
)
def test_vreckon_unit(deg, lat, lon, srange, az, lato, lono):
    if not deg:
        lat, lon, az, lato, lono = map(radians, (lat, lon, az, lato, lono))

    lat1, lon1 = vincenty.vreckon(lat, lon, srange, az, deg=deg)

    assert lat1 == approx(lato)
    assert isinstance(lat1, float)

    assert lon1 == approx(lono, rel=0.001)
    assert isinstance(lon1, float)


def test_negative_lon_stays_negative():
    """Regression test for issue #88: vreckon with negative start longitude
    must return a negative (or near-zero) destination longitude, not 358+."""
    lat2, lon2 = vincenty.vreckon(52.22610277777778, -1.2696583333333333, 839.63, 63.02)
    assert lat2 == approx(52.22952562862266, rel=1e-6)
    assert -180 <= lon2 <= 0, f"lon2={lon2} should be negative for a point west of prime meridian"
    assert lon2 == approx(-1.258707, abs=1e-3)


def test_az_vector():
    np = pytest.importorskip("numpy")
    az = np.array(az1)
    a, b = vincenty.vreckon(*ll0, sr1[0], az)
    assert a == approx(lat2)
    assert b == approx(lon2)


def test_both_vector():
    np = pytest.importorskip("numpy")
    sr = np.array(sr1)
    az = np.array(az1)

    a, b = vincenty.vreckon(10, 20, sr, az)
    assert a == approx(lat3)
    assert b == approx(lon3)
