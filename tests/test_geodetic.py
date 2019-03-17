#!/usr/bin/env python
import pytest
from pytest import approx
from math import radians, nan
import numpy as np

import pymap3d as pm

lla0 = (42, -82, 200)
rlla0 = (radians(lla0[0]), radians(lla0[1]), lla0[2])
lla1 = (42.002582, -81.997752, 1.1397018e3)
rlla1 = (np.radians(lla1[0]), np.radians(lla1[1]), lla1[2])

xyz0 = (660675.2518247,
        -4700948.68316,
        4245737.66222)

aer0 = (33, 70, 1000)
raer0 = (np.radians(aer0[0]), np.radians(aer0[1]), aer0[2])

E = pm.Ellipsoid()

atol_dist = 1e-6  # 1 micrometer


@pytest.mark.parametrize('lla',
                         [(42, -82, 200),
                          ([42], [-82], [200]),
                          (np.array(42), np.array(-82), np.array(200)),
                          (np.array([42]), np.array([-82]), np.array([200])),
                          (np.atleast_3d(42), np.atleast_3d(-82), np.atleast_3d(200))],
                         ids=('scalar', 'list', '0d', '1d', '3d'))
def test_scalar_geodetic2ecef(lla):
    """
    verify we can handle the wide variety of input data type users might use
    """
    x0, y0, z0 = pm.geodetic2ecef(*lla)

    assert (x0, y0, z0) == approx(xyz0)


@pytest.mark.parametrize('xyz',
                         [(xyz0[0], xyz0[1], xyz0[2]),
                          ([xyz0[0]], [xyz0[1]], [xyz0[2]]),
                          (np.array(xyz0[0]), np.array(xyz0[1]), np.array(xyz0[2])),
                          (np.array([xyz0[0]]), np.array([xyz0[1]]), np.array([xyz0[2]])),
                          (np.atleast_3d(xyz0[0]), np.atleast_3d(xyz0[1]), np.atleast_3d(xyz0[2]))],
                         ids=('scalar', 'list', '0d', '1d', '3d'))
def test_scalar_ecef2geodetic(xyz):
    """
    verify we can handle the wide variety of input data type users might use
    """
    lat, lon, alt = pm.ecef2geodetic(*xyz)

    assert [lat, lon, alt] == approx(lla0, rel=1e-4)


@pytest.mark.parametrize('xyz',
                         [(0, E.a, 50),
                          ([0], [E.a], [50]),
                          (np.array(0), np.array(E.a), np.array(50)),
                          (np.array([0]), np.array([E.a]), np.array([50])),
                          (np.atleast_3d(0), np.atleast_3d(E.a), np.atleast_3d(50))],
                         ids=('scalar', 'list', '0d', '1d', '3d'))
def test_scalar_aer_enu(xyz):
    """
    verify we can handle the wide variety of input data type users might use
    """
    enu = pm.ecef2enu(*xyz, 0, 90, -100)

    assert pm.enu2ecef(*enu, 0, 90, -100) == approx([0, E.a, 50])


def test_xarray():
    xarray = pytest.importorskip('xarray')
    xr_lla = xarray.DataArray(list(lla0))

    xyz = pm.geodetic2ecef(*xr_lla)

    assert xyz == approx(xyz0)
    assert isinstance(xyz[0], xarray.DataArray)
# %%
    xr_xyz = xarray.DataArray(list(xyz0))

    lla = pm.ecef2geodetic(*xr_xyz)

    assert lla == approx(lla0)
    assert isinstance(lla[0], float)  # xarrayness is lost, possibly expensive to keep due to isinstance()


def test_pandas():
    pandas = pytest.importorskip('pandas')
    pd_lla = pandas.Series(lla0)

    xyz = pm.geodetic2ecef(*pd_lla)

    assert xyz == approx(xyz0)
    assert isinstance(xyz[0], float)  # series degenerates to scalars by pandas itself
# %% dataframe degenerates to series
    pd_lla = pandas.DataFrame([[*lla0], [*lla0]], columns=['lat', 'lon', 'alt_m'])
    xyz = pm.geodetic2ecef(pd_lla['lat'], pd_lla['lon'], pd_lla['alt_m'])

    assert xyz[0].values == approx(xyz0[0])
    assert xyz[1].values == approx(xyz0[1])
    assert xyz[2].values == approx(xyz0[2])
    assert isinstance(xyz[0], pandas.Series)


def test_ecef():
    xyz = pm.geodetic2ecef(*lla0)

    assert xyz == approx(xyz0)
    assert pm.geodetic2ecef(*rlla0, deg=False) == approx(xyz)

    with pytest.raises(ValueError):
        pm.geodetic2ecef(-100, lla0[1], lla0[2])

    assert pm.ecef2geodetic(*xyz) == approx(lla0)
    assert pm.ecef2geodetic(*xyz, deg=False) == approx(rlla0)

    assert pm.ecef2geodetic((E.a - 1) / np.sqrt(2),
                            (E.a - 1) / np.sqrt(2), 0) == approx([0, 45, -1])


@pytest.mark.parametrize('lla, xyz', [((0, 0, -1), (E.a - 1, 0, 0)),
                                      ((0, 90, -1), (0, E.a - 1, 0)),
                                      ((0, -90, -1), (0, -E.a + 1, 0)),
                                      ((90, 0, -1), (0, 0, E.b - 1)),
                                      ((90, 15, -1), (0, 0, E.b - 1)),
                                      ((-90, 0, -1), (0, 0, -E.b + 1))
                                      ])
def test_geodetic2ecef(lla, xyz):
    assert pm.geodetic2ecef(*lla) == approx(xyz, abs=atol_dist)


@pytest.mark.parametrize('xyz, lla', [((E.a - 1, 0, 0), (0, 0, -1)),
                                      ((0, E.a - 1, 0), (0, 90, -1)),
                                      ((0, 0, E.b - 1), (90, 0, -1)),
                                      ((0, 0, -E.b + 1), (-90, 0, -1)),
                                      ((-E.a + 1, 0, 0), (0, 180, -1)),
                                      ])
def test_ecef2geodetic(xyz, lla):
    assert pm.ecef2geodetic(*xyz) == approx(lla)


def test_aer():
    lla2 = pm.aer2geodetic(*aer0, *lla0)
    rlla2 = pm.aer2geodetic(*raer0, *rlla0, deg=False)

    with pytest.raises(ValueError):
        pm.aer2geodetic(aer0[0], aer0[1], -1, *lla0)

    assert lla2 == approx(lla1)
    assert rlla2 == approx(rlla1)

    assert pm.geodetic2aer(*lla2, *lla0) == approx(aer0)
    assert pm.geodetic2aer(*rlla2, *rlla0, deg=False) == approx(raer0)


def test_allnan():

    anan = np.empty((10, 10))
    anan.fill(nan)
    assert np.isnan(pm.geodetic2aer(anan, anan, anan, *lla0)).all()
    assert np.isnan(pm.aer2geodetic(anan, anan, anan, *lla0)).all()


def test_somenan():
    xyz = np.stack((xyz0, (nan, nan, nan)))

    lat, lon, alt = pm.ecef2geodetic(xyz[:, 0], xyz[:, 1], xyz[:, 2])
    assert (lat[0], lon[0], alt[0]) == approx(lla0)


if __name__ == '__main__':
    pytest.main(['-x', __file__])
