from math import nan

import pymap3d.los as los
import pytest
from pytest import approx

lla0 = (42, -82, 200)


@pytest.mark.parametrize(
    "az,tilt,lat,lon,sr",
    [
        (0, 0, 10, -20, 1e3),
        (0, 90, nan, nan, nan),
        (0, 45, 10.009041667, -20, 1.41432515e3),
        (45, 45, 10.00639336, -19.993549978, 1.414324795e3),
        (125, 45, 9.99481382, -19.992528, 1.414324671e3),
    ],
)
def test_losint(az, tilt, lat, lon, sr):

    lla0 = (10, -20, 1e3)
    lat1, lon1, sr1 = los.lookAtSpheroid(*lla0, az, tilt=tilt)

    nan_ok = True if tilt == 90 else False

    assert lat1 == approx(lat, nan_ok=nan_ok)
    assert lon1 == approx(lon, nan_ok=nan_ok)
    assert sr1 == approx(sr, nan_ok=nan_ok)
    assert isinstance(lat1, float)
    assert isinstance(lon1, float)
    assert isinstance(sr1, float)


def test_badval():

    with pytest.raises(ValueError):
        los.lookAtSpheroid(0, 0, -1, 0, 0)


def test_array_los():
    np = pytest.importorskip("numpy")

    az = [0.0, 10.0, 125.0]
    tilt = [30.0, 45.0, 90.0]

    lat, lon, sr = los.lookAtSpheroid(*lla0, az, tilt)

    truth = np.array(
        [
            [42.00103959, lla0[1], 230.9413173],
            [42.00177328, -81.9995808, 282.84715651],
            [nan, nan, nan],
        ]
    )

    assert np.column_stack((lat, lon, sr)) == approx(truth, nan_ok=True)

    lat, lon, sr = los.lookAtSpheroid([lla0[0]] * 3, [lla0[1]] * 3, [lla0[2]] * 3, az, tilt)
    assert np.column_stack((lat, lon, sr)) == approx(truth, nan_ok=True)


def test_xarray_los():
    xarray = pytest.importorskip("xarray")

    lla = xarray.DataArray(list(lla0))
    az = xarray.DataArray([0.0] * 2)
    tilt = xarray.DataArray([30.0] * 2)

    lat, lon, sr = los.lookAtSpheroid(*lla, az, tilt)
    assert lat == approx(42.00103959)
    assert lon == approx(lla0[1])
    assert sr == approx(230.9413173)


def test_pandas_los():
    pandas = pytest.importorskip("pandas")

    lla = pandas.Series(lla0)
    az = pandas.Series([0.0] * 2)
    tilt = pandas.Series([30.0] * 2)

    lat, lon, sr = los.lookAtSpheroid(*lla, az, tilt)
    assert lat == approx(42.00103959)
    assert lon == approx(lla0[1])
    assert sr == approx(230.9413173)
