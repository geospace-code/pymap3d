import pymap3d.lox as lox
import pytest
from pytest import approx


@pytest.mark.parametrize("lat,dist", [(0, 0), (90, 10001965.729)])
def test_meridian_dist(lat, dist):
    assert lox.meridian_dist(lat) == approx(dist)


@pytest.mark.parametrize(
    "lat1,lat2,arclen",
    [
        (0, 0, 0),
        (0, 90, 10001965.729),
        (0, -90, 10001965.729),
        (0, 40, 4429529.03035058),
        (40, 80, 4455610.84159),
    ],
)
def test_meridian_arc(lat1, lat2, arclen):
    """
    meridianarc(deg2rad(40), deg2rad(80), wgs84Ellipsoid)
    """

    assert lox.meridian_arc(lat1, lat2) == approx(arclen)


@pytest.mark.parametrize(
    "lon1,lon2,lat,dist",
    [
        (0, 0, 0, 0),
        (0, 90, 0, 10018754.1714),
        (0, -90, 0, 10018754.1714),
        (90, 0, 0, 10018754.1714),
        (-90, 0, 0, 10018754.1714),
    ],
)
def test_departure(lon1, lon2, lat, dist):
    assert lox.departure(lon1, lon2, lat) == approx(dist)


@pytest.mark.parametrize(
    "lat1,lon1,lat2,lon2,arclen,az",
    [
        (40, -80, 65, -148, 5248666.20853187, 302.0056736),
        (0, 0, 0, 1, 111319.49, 90),
        (0, 0, 0, -1, 111319.49, 270),
        (0, 1, 0, 0, 111319.49, 270),
        (0, -1, 0, 0, 111319.49, 90),
        (1, 0, 0, 0, 110574.4, 180),
        (-1, 0, 0, 0, 110574.4, 0),
    ],
)
def test_loxodrome_inverse(lat1, lon1, lat2, lon2, arclen, az):
    """
    distance('rh', 40, -80, 65, -148, wgs84Ellipsoid)
    azimuth('rh', 40, -80, 65, -148, wgs84Ellipsoid)
    """
    rhdist, rhaz = lox.loxodrome_inverse(lat1, lon1, lat2, lon2)

    assert rhdist == approx(arclen)
    assert rhaz == approx(az)
    assert isinstance(rhdist, float)
    assert isinstance(rhaz, float)


def test_numpy_loxodrome_inverse():
    pytest.importorskip("numpy")
    d, a = lox.loxodrome_inverse([40, 40], [-80, -80], 65, -148)
    assert d == approx(5248666.209)
    assert a == approx(302.00567)

    d, a = lox.loxodrome_inverse([40, 40], [-80, -80], [65, 65], -148)
    d, a = lox.loxodrome_inverse([40, 40], [-80, -80], 65, [-148, -148])


def test_numpy_2d_loxodrome_inverse():
    pytest.importorskip("numpy")
    d, a = lox.loxodrome_inverse([[40, 40], [40, 40]], [[-80, -80], [-80, -80]], 65, -148)
    assert d == approx(5248666.209)
    assert a == approx(302.00567)

    d, a = lox.loxodrome_inverse(
        [[40, 40], [40, 40]], [[-80, -80], [-80, -80]], [[65, 65], [65, 65]], -148
    )
    d, a = lox.loxodrome_inverse(
        [[40, 40], [40, 40]], [[-80, -80], [-80, -80]], 65, [[-148, -148], [-148, -148]]
    )
    d, a = lox.loxodrome_inverse(40, -80, [[65, 65], [65, 65]], [[-148, -148], [-148, -148]])


@pytest.mark.parametrize(
    "lat0,lon0,rng,az,lat1,lon1",
    [
        (40, -80, 10000, 30, 40.0000779959676, -79.9999414477481),
        (35, 140, 50000, 90, 35, 140.548934481815),
        (35, 140, 50000, -270, 35, 140.548934481815),
        (35, 140, 50000, 270, 35, 139.451065518185),
        (35, 140, 50000, -90, 35, 139.451065518185),
        (0, 0, 0, 0, 0, 0),
        (0, 0, 10018754.17, 90, 0, 90),
        (0, 0, 10018754.17, -90, 0, -90),
        (0, 0, 110574.4, 180, -1, 0),
        (-1, 0, 110574.4, 0, 0, 0),
    ],
)
def test_loxodrome_direct(lat0, lon0, rng, az, lat1, lon1):
    """
    reckon('rh', 40, -80, 10, 30, wgs84Ellipsoid)
    """
    lat2, lon2 = lox.loxodrome_direct(lat0, lon0, rng, az)
    assert lat2 == approx(lat1, rel=0.005, abs=1e-6)
    assert lon2 == approx(lon1, rel=0.001)
    assert isinstance(lat2, float)
    assert isinstance(lon2, float)


def test_numpy_loxodrome_direct():
    pytest.importorskip("numpy")
    lat, lon = lox.loxodrome_direct([40, 40], [-80, -80], [10, 10], [30, 30])
    assert lat == approx(40.000078)
    assert lon == approx(-79.99994145)

    lat, lon = lox.loxodrome_direct([40, 40], [-80, -80], 10, 30)
    lat, lon = lox.loxodrome_direct([40, 40], [-80, -80], [10, 10], 30)
    lat, lon = lox.loxodrome_direct([40, 40], [-80, -80], 10, [30, 30])


@pytest.mark.parametrize("lat,lon", [([0, 45, 90], [0, 45, 90])])
def test_meanm(lat, lon):
    pytest.importorskip("numpy")
    assert lox.meanm(lat, lon) == approx([47.26967, 18.460557])
