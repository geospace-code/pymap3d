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
def test_unit_vdist(lat, lon, lat1, lon1, srange, az):
    dist, az1 = vincenty.vdist(lat, lon, lat1, lon1)
    assert dist == approx(srange, rel=0.005)
    assert az1 == approx(az)

    assert isinstance(dist, float)
    assert isinstance(az1, float)


def test_vector():
    np = pytest.importorskip("numpy")
    lat2 = np.array([10.02137267, 10.01917819])
    lon2 = np.array([20.0168471, 20.0193493])
    asr, aaz = vincenty.vdist(10, 20, lat2, lon2)

    assert 3e3 == approx(asr)
    assert aaz == approx([38, 45])


@pytest.mark.parametrize("shape", [(3,), (1, 3), (3, 1)])
@pytest.mark.parametrize("deg", [True, False])
@pytest.mark.parametrize("error_mode", ["warn", "raise", "ignore"])
def test_vector_coincident_pair(shape, deg, error_mode):
    np = pytest.importorskip("numpy")
    coords = np.array(
        [[0, 10, 40], [0, 20, -70], [0, 11, 41], [0, 21, -72]], dtype=float
    )
    if not deg:
        coords = np.radians(coords)
    expected = np.array([vincenty.vdist(*point, deg=deg) for point in coords.T]).T
    original = coords.copy()

    with np.errstate(all=error_mode):
        policy = np.geterr().copy()
        distance, azimuth = vincenty.vdist(
            *(value.reshape(shape) for value in coords), deg=deg
        )
        assert np.geterr() == policy

    assert distance.shape == azimuth.shape == shape
    np.testing.assert_allclose(distance.ravel(), expected[0], rtol=0, atol=1e-5)
    np.testing.assert_allclose(azimuth.ravel(), expected[1], rtol=0, atol=1e-8)
    np.testing.assert_array_equal(coords, original)


def test_vector_coincident_scalar_origin():
    np = pytest.importorskip("numpy")
    distance, azimuth = vincenty.vdist(
        0.0, 0.0, np.array([0.0, 10.0]), np.array([0.0, 20.0])
    )
    expected = vincenty.vdist(0.0, 0.0, 10.0, 20.0)

    assert distance[0] == azimuth[0] == 0.0
    assert distance[1] == approx(expected[0], rel=1e-12)
    assert azimuth[1] == approx(expected[1], abs=1e-10)


@pytest.mark.parametrize("lat,lon,slantrange,az", [(10, 20, 3e3, 38), (0, 0, 0, 0)])
def test_identity(lat, lon, slantrange, az):
    lat1, lon1 = vincenty.vreckon(lat, lon, slantrange, az)

    dist, az1 = vincenty.vdist(lat, lon, lat1, lon1)

    assert dist == approx(slantrange)
    assert az1 == approx(az)
