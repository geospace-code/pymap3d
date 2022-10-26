import pymap3d.vincenty as vincenty
import pytest
from pytest import approx


@pytest.mark.parametrize(
    "lat,lon,lat1,lon1,srange,az",
    [
        (0, 0, 0, 0, 0, 0),
        (0, 0, 0, 90, 1.001875e7, 90),
        (0, 0, 0, -90, 1.001875e7, 270),
        (0, 0, 0, 180, 2.00375e7, 90),
        (0, 0, 0, -180, 2.00375e7, 90),
        (0, 0, 0, 4, 445277.96, 90),
        (0, 0, 0, 5, 556597.45, 90),
        (0, 0, 0, 6, 667916.94, 90),
        (0, 0, 0, -6, 667916.94, 270),
        (0, 0, 0, 7, 779236.44, 90),
        (1e-16, 1e-16, 1e-16, 1, 111319.49, 90),
        (90, 0, 0, 0, 1.00019657e7, 180),
        (90, 0, -90, 0, 2.000393145e7, 180),
    ],
)
def test_unit(lat, lon, lat1, lon1, srange, az):
    dist, az1 = vincenty.vdist(lat, lon, lat1, lon1)
    assert dist == approx(srange, rel=0.005)
    assert az1 == approx(az)

    assert isinstance(dist, float)
    assert isinstance(az1, float)


def test_vector():
    pytest.importorskip("numpy")
    asr, aaz = vincenty.vdist(10, 20, [10.02137267, 10.01917819], [20.0168471, 20.0193493])

    assert 3e3 == approx(asr)
    assert aaz == approx([38, 45])


@pytest.mark.parametrize("lat,lon,slantrange,az", [(10, 20, 3e3, 38), (0, 0, 0, 0)])
def test_identity(lat, lon, slantrange, az):
    lat1, lon1 = vincenty.vreckon(lat, lon, slantrange, az)

    dist, az1 = vincenty.vdist(lat, lon, lat1, lon1)

    assert dist == approx(slantrange)
    assert az1 == approx(az)
