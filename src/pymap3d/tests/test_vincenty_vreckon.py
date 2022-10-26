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


@pytest.mark.parametrize(
    "lat,lon,srange,az,lato,lono",
    [
        (0, 0, 0, 0, 0, 0),
        (0, 0, 1.001875e7, 90, 0, 90),
        (0, 0, 1.001875e7, 270, 0, 270),
        (0, 0, 1.001875e7, -90, 0, 270),
        (0, 0, 2.00375e7, 90, 0, 180),
        (0, 0, 2.00375e7, 270, 0, 180),
        (0, 0, 2.00375e7, -90, 0, 180),
    ],
)
def test_unit(lat, lon, srange, az, lato, lono):
    lat1, lon1 = vincenty.vreckon(lat, lon, srange, az)

    assert lat1 == approx(lato)
    assert isinstance(lat1, float)

    assert lon1 == approx(lono, rel=0.001)
    assert isinstance(lon1, float)


def test_az_vector():
    pytest.importorskip("numpy")
    a, b = vincenty.vreckon(*ll0, sr1[0], az1)
    assert a == approx(lat2)
    assert b == approx(lon2)


def test_both_vector():
    pytest.importorskip("numpy")
    a, b = vincenty.vreckon(10, 20, sr1, az1)
    assert a == approx(lat3)
    assert b == approx(lon3)
