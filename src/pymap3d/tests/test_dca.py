import pytest

import pymap3d as pm


def get_ellipsoid_params():
    ell = pm.Ellipsoid.from_name("wgs84")
    return ell.semimajor_axis, ell.semiminor_axis


A, B = get_ellipsoid_params()


@pytest.mark.parametrize(
    "enu,dca,heading",
    [
        ((0, 0, 0), (0, 0, 0), 15),
        ((-7.0710678118654755, 12.24744871391589, 1000), (10, 10, 1000), 15),
        ((-2.455756079379457, 13.927284806400378, 1000), (10, 10, 1000), 35),
    ],
)
def test_enu_dca(enu, dca, heading):

    assert pm.dca2enu(*dca, heading=heading) == pytest.approx(enu)
    assert pm.enu2dca(*enu, heading=heading) == pytest.approx(dca)

    ned = enu[1], enu[0], -enu[2]
    assert pm.dca2ned(*dca, heading=heading) == pytest.approx(ned)
    assert pm.ned2dca(*ned, heading=heading) == pytest.approx(dca)


def test_enu_dca_numpy():
    np = pytest.importorskip("numpy")
    e = np.array([0, 0, 0])
    n = np.array([0, 0, 0])
    u = np.array([0, 0, 0])
    d = np.array([0, 0, 0])
    c = np.array([0, 0, 0])
    a = np.array([0, 0, 0])
    heading = 15
    assert pm.dca2enu(d, c, a, heading=heading) == pytest.approx((0, 0, 0))
    assert pm.enu2dca(e, n, u, heading=heading) == pytest.approx((0, 0, 0))


@pytest.mark.parametrize(
    "ecef,dca,heading",
    [
        ((6027079.293014112, 1614951.0292814733, 1317408.7685803245), (0, 0, 0), 15),
        (
            (6028023.481548891, 1615196.7033287943, 1317628.660083717),
            (10, 10, 1000),
            15,
        ),
    ],
)
def test_ecef_dca(ecef, dca, heading):
    lat0, lon0, h0 = 12.0, 15.0, 30.0

    assert pm.dca2ecef(*dca, lat0, lon0, h0, heading) == pytest.approx(ecef)
    assert pm.ecef2dca(*ecef, lat0, lon0, h0, heading) == pytest.approx(dca, abs=1e-9)


def test_geodetic_dca():
    lat0, lon0, h0 = 12.0, 15.0, 30.0
    heading = 15.0

    lat, lon, h = 12.1, 15.1, 30.1
    dca = pm.geodetic2dca(lat, lon, h, lat0, lon0, h0, heading)

    assert pm.dca2geodetic(*dca, lat0, lon0, h0, heading) == pytest.approx((lat, lon, h))


def test_aer_dca():
    heading = 15.0

    az, el, r = 10.0, 20.0, 1000.0
    dca = pm.aer2dca(az, el, r, heading)

    assert pm.dca2aer(*dca, heading) == pytest.approx((az, el, r))
