#!/usr/bin/env python
import pytest
from pytest import approx
from math import radians, nan, sqrt, isnan

import pymap3d as pm

lla0 = (42, -82, 200)
rlla0 = (radians(lla0[0]), radians(lla0[1]), lla0[2])
lla1 = (42.002582, -81.997752, 1.1397018e3)
rlla1 = (radians(lla1[0]), radians(lla1[1]), lla1[2])

xyz0 = (660675.2518247,
        -4700948.68316,
        4245737.66222)

aer0 = (33, 70, 1000)
raer0 = (radians(aer0[0]), radians(aer0[1]), aer0[2])

ELL = pm.Ellipsoid()
A = ELL.semimajor_axis
B = ELL.semiminor_axis

atol_dist = 1e-6  # 1 micrometer


@pytest.mark.parametrize('lla',
                         [(42, -82, 200),
                          ([42], [-82], [200])],
                         ids=('scalar', 'list'))
def test_scalar_geodetic2ecef(lla):
    """
    verify we can handle the wide variety of input data type users might use
    """
    x0, y0, z0 = pm.geodetic2ecef(*lla)

    assert (x0, y0, z0) == approx(xyz0)


def test_3d_geodetic2ecef():
    np = pytest.importorskip("numpy")
    lla = (np.atleast_3d(42), np.atleast_3d(-82), np.atleast_3d(200))
    x0, y0, z0 = pm.geodetic2ecef(*lla)

    assert (x0, y0, z0) == approx(xyz0)


@pytest.mark.parametrize('xyz',
                         [(xyz0[0], xyz0[1], xyz0[2]),
                          ([xyz0[0]], [xyz0[1]], [xyz0[2]])],
                         ids=('scalar', 'list'))
def test_scalar_ecef2geodetic(xyz):
    """
    verify we can handle the wide variety of input data type users might use
    """
    lat, lon, alt = pm.ecef2geodetic(*xyz)

    assert [lat, lon, alt] == approx(lla0, rel=1e-4)


def test_3d_ecef2geodetic():
    np = pytest.importorskip("numpy")
    xyz = (np.atleast_3d(xyz0[0]), np.atleast_3d(xyz0[1]), np.atleast_3d(xyz0[2]))

    lat, lon, alt = pm.ecef2geodetic(*xyz)

    assert [lat, lon, alt] == approx(lla0, rel=1e-4)


@pytest.mark.parametrize('xyz',
                         [(0, A, 50),
                          ([0], [A], [50])],
                         ids=('scalar', 'list'))
def test_scalar_aer_enu(xyz):
    """
    verify we can handle the wide variety of input data type users might use
    """
    enu = pm.ecef2enu(*xyz, 0, 90, -100)

    assert pm.enu2ecef(*enu, 0, 90, -100) == approx([0, A, 50])


def test_3d_aer_enu():
    np = pytest.importorskip("numpy")
    xyz = (np.atleast_3d(0), np.atleast_3d(A), np.atleast_3d(50))

    enu = pm.ecef2enu(*xyz, 0, 90, -100)
    assert pm.enu2ecef(*enu, 0, 90, -100) == approx([0, A, 50])


def test_xarray():
    xarray = pytest.importorskip('xarray')
    xr_lla = xarray.DataArray(list(lla0))

    xyz = pm.geodetic2ecef(*xr_lla)

    assert xyz == approx(xyz0)
# %%
    xr_xyz = xarray.DataArray(list(xyz0))

    lla = pm.ecef2geodetic(*xr_xyz)

    assert lla == approx(lla0)


def test_pandas():
    pandas = pytest.importorskip('pandas')
    pd_lla = pandas.Series(lla0)

    xyz = pm.geodetic2ecef(*pd_lla)

    assert xyz == approx(xyz0)
# %% dataframe degenerates to series
    pd_lla = pandas.DataFrame([[*lla0], [*lla0]], columns=['lat', 'lon', 'alt_m'])
    xyz = pm.geodetic2ecef(pd_lla['lat'], pd_lla['lon'], pd_lla['alt_m'])

    assert xyz[0] == approx(xyz0[0])
    assert xyz[1] == approx(xyz0[1])
    assert xyz[2] == approx(xyz0[2])


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

    assert pm.ecef2geodetic((A - 1) / sqrt(2),
                            (A - 1) / sqrt(2), 0) == approx([0, 45, -1])


@pytest.mark.parametrize('lla, xyz', [((0, 0, -1), (A - 1, 0, 0)),
                                      ((0, 90, -1), (0, A - 1, 0)),
                                      ((0, -90, -1), (0, -A + 1, 0)),
                                      ((90, 0, -1), (0, 0, B - 1)),
                                      ((90, 15, -1), (0, 0, B - 1)),
                                      ((-90, 0, -1), (0, 0, -B + 1))
                                      ])
def test_geodetic2ecef(lla, xyz):
    assert pm.geodetic2ecef(*lla) == approx(xyz, abs=atol_dist)


@pytest.mark.parametrize('xyz, lla', [((A - 1, 0, 0), (0, 0, -1)),
                                      ((0, A - 1, 0), (0, 90, -1)),
                                      ((0, 0, B - 1), (90, 0, -1)),
                                      ((0, 0, -B + 1), (-90, 0, -1)),
                                      ((-A + 1, 0, 0), (0, 180, -1)),
                                      ])
def test_ecef2geodetic(xyz, lla):
    lat, lon, alt = pm.ecef2geodetic(*xyz)
    assert lat == approx(lla[0])
    assert lon == approx(lla[1])
    assert alt == approx(lla[2])


def test_aer():
    lla2 = pm.aer2geodetic(*aer0, *lla0)
    rlla2 = pm.aer2geodetic(*raer0, *rlla0, deg=False)

    with pytest.raises(ValueError):
        pm.aer2geodetic(aer0[0], aer0[1], -1, *lla0)

    assert lla2 == approx(lla1)
    assert rlla2 == approx(rlla1)

    assert pm.geodetic2aer(*lla2, *lla0) == approx(aer0)
    assert pm.geodetic2aer(*rlla2, *rlla0, deg=False) == approx(raer0)


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


if __name__ == '__main__':
    pytest.main(['-x', __file__])
