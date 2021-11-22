import pytest
from pytest import approx
from math import radians, nan, sqrt, isnan

import pymap3d as pm

lla0 = (42, -82, 200)
rlla0 = (radians(lla0[0]), radians(lla0[1]), lla0[2])

xyz0 = (660675.2518247, -4700948.68316, 4245737.66222)

ELL = pm.Ellipsoid()
A = ELL.semimajor_axis
B = ELL.semiminor_axis

xyzlla = [
    ((A - 1, 0, 0), (0, 0, -1)),
    ((0, A - 1, 0), (0, 90, -1)),
    ((0, 0, B - 1), (90, 0, -1)),
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


@pytest.mark.parametrize("lla", [(42, -82, 200), ([42], [-82], [200])], ids=("scalar", "list"))
def test_scalar_geodetic2ecef(lla):
    """
    verify we can handle the wide variety of input data type users might use
    """
    if isinstance(lla[0], list):
        pytest.importorskip("numpy")

    x0, y0, z0 = pm.geodetic2ecef(*lla)

    assert (x0, y0, z0) == approx(xyz0)


def test_3d_geodetic2ecef():
    np = pytest.importorskip("numpy")
    lla = (np.atleast_3d(42), np.atleast_3d(-82), np.atleast_3d(200))
    x0, y0, z0 = pm.geodetic2ecef(*lla)

    assert (x0, y0, z0) == approx(xyz0)


def test_scalar_ecef2geodetic():
    """
    verify we can handle the wide variety of input data type users might use
    """
    lat, lon, alt = pm.ecef2geodetic(xyz0[0], xyz0[1], xyz0[2])

    assert [lat, lon, alt] == approx(lla0, rel=1e-4)


def test_3d_ecef2geodetic():
    np = pytest.importorskip("numpy")
    xyz = (np.atleast_3d(xyz0[0]), np.atleast_3d(xyz0[1]), np.atleast_3d(xyz0[2]))

    lat, lon, alt = pm.ecef2geodetic(*xyz)

    assert [lat, lon, alt] == approx(lla0, rel=1e-4)


def test_array_ecef2geodetic():
    """
    tests ecef2geodetic can handle numpy array data in addition to singular floats
    """
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


def test_xarray():
    xarray = pytest.importorskip("xarray")
    xr_lla = xarray.DataArray(list(lla0))

    xyz = pm.geodetic2ecef(*xr_lla)

    assert xyz == approx(xyz0)

    xr_xyz = xarray.DataArray(list(xyz0))

    lla = pm.ecef2geodetic(*xr_xyz)

    assert lla == approx(lla0)


def test_pandas():
    pandas = pytest.importorskip("pandas")
    pd_lla = pandas.Series(lla0)

    xyz = pm.geodetic2ecef(*pd_lla)

    assert xyz == approx(xyz0)
    # %% dataframe degenerates to series
    pd_lla = pandas.DataFrame([[*lla0], [*lla0]], columns=["lat", "lon", "alt_m"])
    xyz = pm.geodetic2ecef(pd_lla["lat"], pd_lla["lon"], pd_lla["alt_m"])

    assert xyz[0].values == approx(xyz0[0])
    assert xyz[1].values == approx(xyz0[1])
    assert xyz[2].values == approx(xyz0[2])


def test_ecef():
    xyz = pm.geodetic2ecef(*lla0)

    assert xyz == approx(xyz0)
    x, y, z = pm.geodetic2ecef(*rlla0, deg=False)
    assert x == approx(xyz[0])
    assert y == approx(xyz[1])
    assert z == approx(xyz[2])

    with pytest.raises(ValueError):
        pm.geodetic2ecef(-100, lla0[1], lla0[2])

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
    assert alt == approx(lla[2])


@pytest.mark.parametrize(
    "aer,lla,lla0",
    [
        ((33, 77, 1000), (42.0016981935, -81.99852, 1174.374035), (42, -82, 200)),
        ((0, 90, 10000), (0, 0, 10000), (0, 0, 0)),
    ],
)
def test_aer_geodetic(aer, lla, lla0):
    lat1, lon1, alt1 = pm.aer2geodetic(*aer, *lla0)
    assert lat1 == approx(lla[0])
    assert lon1 == approx(lla[1])
    assert alt1 == approx(lla[2])
    assert isinstance(lat1, float)
    assert isinstance(lon1, float)
    assert isinstance(alt1, float)

    raer = (radians(aer[0]), radians(aer[1]), aer[2])
    rlla0 = (radians(lla0[0]), radians(lla0[1]), lla0[2])
    assert pm.aer2geodetic(*raer, *rlla0, deg=False) == approx(
        (radians(lla[0]), radians(lla[1]), lla[2])
    )

    with pytest.raises(ValueError):
        pm.aer2geodetic(aer[0], aer[1], -1, *lla0)

    assert pm.geodetic2aer(*lla, *lla0) == approx(aer, rel=1e-3)
    assert pm.geodetic2aer(radians(lla[0]), radians(lla[1]), lla[2], *rlla0, deg=False) == approx(
        raer, rel=1e-3
    )


def test_scalar_nan():
    a, e, r = pm.geodetic2aer(nan, nan, nan, *lla0)
    assert isnan(a) and isnan(e) and isnan(r)

    lat, lon, alt = pm.aer2geodetic(nan, nan, nan, *lla0)
    assert isnan(lat) and isnan(lon) and isnan(alt)


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
    lat, lon, alt = pm.ecef2geodetic(
        *np.array(
            [
                [xyz],
            ],
            dtype=np.float32,
        ).T
    )
    assert lat[0] == approx(lla[0])
    assert lon[0] == approx(lla[1])
    assert alt[0] == approx(lla[2])
