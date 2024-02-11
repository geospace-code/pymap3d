from math import isnan, nan, radians, sqrt

import pymap3d as pm
import pytest
from pytest import approx

lla0 = (42, -82, 200)
rlla0 = (radians(lla0[0]), radians(lla0[1]), lla0[2])

xyz0 = (660675.2518247, -4700948.68316, 4245737.66222)

ELL = pm.Ellipsoid.from_name("wgs84")
A = ELL.semimajor_axis
B = ELL.semiminor_axis

xyzlla = [
    ((A, 0, 0), (0, 0, 0)),
    ((A - 1, 0, 0), (0, 0, -1)),
    ((A + 1, 0, 0), (0, 0, 1)),
    ((0.1 * A, 0, 0), (0, 0, -0.9 * A)),
    ((0.001 * A, 0, 0), (0, 0, -0.999 * A)),
    ((0, A, 0), (0, 90, 0)),
    ((0, A - 1, 0), (0, 90, -1)),
    ((0, A + 1, 0), (0, 90, 1)),
    ((0, 0.1 * A, 0), (0, 90, -0.9 * A)),
    ((0, 0.001 * A, 0), (0, 90, -0.999 * A)),
    ((0, 0, B), (90, 0, 0)),
    ((0, 0, B + 1), (90, 0, 1)),
    ((0, 0, B - 1), (90, 0, -1)),
    ((0, 0, 0.1 * B), (90, 0, -0.9 * B)),
    ((0, 0, 0.001 * B), (90, 0, -0.999 * B)),
    ((0, 0, B - 1), (89.999999, 0, -1)),
    ((0, 0, B - 1), (89.99999, 0, -1)),
    ((0, 0, -B + 1), (-90, 0, -1)),
    ((0, 0, -B + 1), (-89.999999, 0, -1)),
    ((0, 0, -B + 1), (-89.99999, 0, -1)),
    ((-A + 1, 0, 0), (0, 180, -1)),
]

llaxyz = [
    ((0, 0, -1), (A - 1, 0, 0)),
    ((0, 90, -1), (0, A - 1, 0)),
    ((0, -90, -1), (0, -A + 1, 0)),
    ((90, 0, -1), (0, 0, B - 1)),
    ((90, 15, -1), (0, 0, B - 1)),
    ((-90, 0, -1), (0, 0, -B + 1)),
]


atol_dist = 1e-6  # 1 micrometer


@pytest.mark.parametrize("lla", [lla0, ([lla0[0]], [lla0[1]], [lla0[2]])], ids=("scalar", "list"))
def test_scalar_geodetic2ecef(lla):
    """
    verify we can handle the wide variety of input data type users might use
    """

    if isinstance(lla[0], list):
        np = pytest.importorskip("numpy")
        scalar = False
    else:
        scalar = True

    xyz = pm.geodetic2ecef(*lla)
    lla1 = pm.ecef2geodetic(*xyz)

    try:
        np.testing.assert_allclose(lla1, lla, rtol=1e-4)
    except NameError:
        assert lla1 == approx(lla, rel=1e-4)

    if scalar:
        assert all(isinstance(n, float) for n in xyz)
        assert all(isinstance(n, float) for n in lla1)


def test_array_geodetic2ecef():
    np = pytest.importorskip("numpy")

    lla = (np.asarray(lla0[0]), np.asarray(lla0[1]), np.asarray(lla0[2]))
    xyz = pm.geodetic2ecef(*lla)
    np.testing.assert_allclose(pm.ecef2geodetic(*xyz), lla)

    lla = (np.atleast_1d(lla0[0]), np.atleast_1d(lla0[1]), np.atleast_1d(lla0[2]))
    xyz = pm.geodetic2ecef(*lla)
    np.testing.assert_allclose(pm.ecef2geodetic(*xyz), lla)


@pytest.mark.parametrize("xyz", [xyz0, ([xyz0[0]], [xyz0[1]], [xyz0[2]])], ids=("scalar", "list"))
def test_scalar_ecef2geodetic(xyz):
    """
    verify we can handle the wide variety of input data type users might use
    """

    if isinstance(xyz[0], list):
        np = pytest.importorskip("numpy")
        scalar = False
    else:
        scalar = True

    lla = pm.ecef2geodetic(*xyz)
    xyz1 = pm.geodetic2ecef(*lla)

    try:
        np.testing.assert_allclose(xyz1, xyz, rtol=1e-4)
    except NameError:
        assert xyz1 == approx(xyz, rel=1e-4)

    if scalar:
        assert all(isinstance(n, float) for n in xyz1)
        assert all(isinstance(n, float) for n in lla)


def test_array_ecef2geodetic():
    np = pytest.importorskip("numpy")

    xyz = (np.asarray(xyz0[0]), np.asarray(xyz0[1]), np.asarray(xyz0[2]))
    lla = pm.ecef2geodetic(*xyz)
    np.testing.assert_allclose(pm.geodetic2ecef(*lla), xyz)

    xyz = (np.atleast_1d(xyz0[0]), np.atleast_1d(xyz0[1]), np.atleast_1d(xyz0[2]))
    lla = pm.ecef2geodetic(*xyz)
    np.testing.assert_allclose(pm.geodetic2ecef(*lla), xyz)


def test_inside_ecef2geodetic():
    np = pytest.importorskip("numpy")
    # test values with no points inside ellipsoid
    lla0_array = (
        np.array([lla0[0], lla0[0]]),
        np.array([lla0[1], lla0[1]]),
        np.array([lla0[2], lla0[2]]),
    )
    xyz = pm.geodetic2ecef(*lla0_array)
    lats, lons, alts = pm.ecef2geodetic(*xyz)

    assert lats == approx(lla0_array[0])
    assert lons == approx(lla0_array[1])
    assert alts == approx(lla0_array[2])

    # test values with some (but not all) points inside ellipsoid
    lla0_array_inside = (
        np.array([lla0[0], lla0[0]]),
        np.array([lla0[1], lla0[1]]),
        np.array([lla0[2], -lla0[2]]),
    )
    xyz = pm.geodetic2ecef(*lla0_array_inside)
    lats, lons, alts = pm.ecef2geodetic(*xyz)

    assert lats == approx(lla0_array_inside[0])
    assert lons == approx(lla0_array_inside[1])
    assert alts == approx(lla0_array_inside[2])


def test_xarray_ecef():
    xarray = pytest.importorskip("xarray")

    lla = xarray.DataArray(list(lla0))

    xyz = pm.geodetic2ecef(*lla)
    lla1 = pm.ecef2geodetic(*xyz)
    assert lla1 == approx(lla)


def test_pandas_ecef():
    pandas = pytest.importorskip("pandas")

    x, y, z = pm.geodetic2ecef(
        pandas.Series(lla0[0]), pandas.Series(lla0[1]), pandas.Series(lla0[2])
    )

    lat, lon, alt = pm.ecef2geodetic(pandas.Series(x), pandas.Series(y), pandas.Series(z))
    assert lat == approx(lla0[0])
    assert lon == approx(lla0[1])
    assert alt == approx(lla0[2])


def test_ecef():
    xyz = pm.geodetic2ecef(*lla0)

    assert xyz == approx(xyz0)
    x, y, z = pm.geodetic2ecef(*rlla0, deg=False)
    assert x == approx(xyz[0])
    assert y == approx(xyz[1])
    assert z == approx(xyz[2])

    assert pm.ecef2geodetic(*xyz) == approx(lla0)
    assert pm.ecef2geodetic(*xyz, deg=False) == approx(rlla0)

    assert pm.ecef2geodetic((A - 1) / sqrt(2), (A - 1) / sqrt(2), 0) == approx([0, 45, -1])


@pytest.mark.parametrize("lla, xyz", llaxyz)
def test_geodetic2ecef(lla, xyz):
    assert pm.geodetic2ecef(*lla) == approx(xyz, abs=atol_dist)


@pytest.mark.parametrize("xyz, lla", xyzlla)
def test_ecef2geodetic(xyz, lla):
    lat, lon, alt = pm.ecef2geodetic(*xyz)
    assert lat == approx(lla[0])
    assert lon == approx(lla[1])
    assert alt == approx(lla[2], abs=1e-9)


@pytest.mark.parametrize(
    "aer,lla,lla0",
    [
        ((33, 77, 1000), (42.0016981935, -81.99852, 1174.374035), (42, -82, 200)),
        ((0, 90, 10000), (0, 0, 10000), (0, 0, 0)),
    ],
)
def test_aer_geodetic(aer, lla, lla0):
    lla1 = pm.aer2geodetic(*aer, *lla0)
    assert lla1 == approx(lla)
    assert all(isinstance(n, float) for n in lla1)

    raer = (radians(aer[0]), radians(aer[1]), aer[2])
    rlla0 = (radians(lla0[0]), radians(lla0[1]), lla0[2])
    lla1 = pm.aer2geodetic(*raer, *rlla0, deg=False)
    assert lla1 == approx((radians(lla[0]), radians(lla[1]), lla[2]))
    assert all(isinstance(n, float) for n in lla1)

    with pytest.raises(ValueError):
        pm.aer2geodetic(aer[0], aer[1], -1, *lla0)

    assert pm.geodetic2aer(*lla, *lla0) == approx(aer, rel=1e-3)
    assert pm.geodetic2aer(radians(lla[0]), radians(lla[1]), lla[2], *rlla0, deg=False) == approx(
        raer, rel=1e-3
    )


def test_scalar_nan():
    aer = pm.geodetic2aer(nan, nan, nan, *lla0)
    assert all(isnan(n) for n in aer)

    llat = pm.aer2geodetic(nan, nan, nan, *lla0)
    assert all(isnan(n) for n in llat)


def test_allnan():
    np = pytest.importorskip("numpy")
    anan = np.empty((10, 10))
    anan.fill(nan)
    assert np.isnan(pm.geodetic2aer(anan, anan, anan, *lla0)).all()
    assert np.isnan(pm.aer2geodetic(anan, anan, anan, *lla0)).all()


def test_somenan():
    np = pytest.importorskip("numpy")
    xyz = np.stack((xyz0, (nan, nan, nan)))

    lat, lon, alt = pm.ecef2geodetic(xyz[:, 0], xyz[:, 1], xyz[:, 2])
    assert (lat[0], lon[0], alt[0]) == approx(lla0)


@pytest.mark.parametrize("xyz, lla", xyzlla)
def test_numpy_ecef2geodetic(xyz, lla):
    np = pytest.importorskip("numpy")
    lla1 = pm.ecef2geodetic(
        *np.array(
            [
                [xyz],
            ],
            dtype=np.float32,
        ).T
    )
    assert lla1 == approx(lla)


@pytest.mark.parametrize("lla, xyz", llaxyz)
def test_numpy_geodetic2ecef(lla, xyz):
    np = pytest.importorskip("numpy")
    xyz1 = pm.geodetic2ecef(
        *np.array(
            [
                [lla],
            ],
            dtype=np.float32,
        ).T
    )

    atol_dist = 1  # meters
    assert xyz1 == approx(xyz, abs=atol_dist)
