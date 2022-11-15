import pymap3d as pm
import pytest
from pymap3d.vincenty import vreckon
from pytest import approx

lla0 = [42, -82, 200]


def test_compare_vicenty() -> None:
    taz, tsr = 38, 3000
    pyproj = pytest.importorskip("pyproj")

    lat2, lon2 = vreckon(10, 20, tsr, taz)

    p4lon, p4lat, p4a21 = pyproj.Geod(ellps="WGS84").fwd(lon2, lat2, taz, tsr)
    assert p4lon == approx(lon2, rel=0.0025)
    assert p4lat == approx(lat2, rel=0.0025)

    p4az, p4a21, p4sr = pyproj.Geod(ellps="WGS84").inv(20, 10, lon2, lat2)
    assert (p4az, p4sr) == approx((taz, tsr))


def test_compare_geodetic() -> None:
    pyproj = pytest.importorskip("pyproj")

    xyz = pm.geodetic2ecef(*lla0)  # type: ignore[call-overload]

    ecef = pyproj.Proj(proj="geocent", ellps="WGS84", datum="WGS84")
    lla = pyproj.Proj(proj="latlong", ellps="WGS84", datum="WGS84")

    assert pyproj.transform(lla, ecef, lla0[1], lla0[0], lla0[2]) == approx(xyz)
    assert pyproj.transform(ecef, lla, *xyz) == approx((lla0[1], lla0[0], lla0[2]))
