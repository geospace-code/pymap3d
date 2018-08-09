import pytest
from pymap3d.vincenty import vreckon, vdist
from numpy.testing import assert_allclose
import pymap3d as pm

az = 38
sr = 3e3
lla0 = (42, -82, 200)


def test_vincenty():

    lat2, lon2, a21 = vreckon(10, 20, sr, az)
    assert_allclose((lat2, lon2, a21),
                    (10.02137267, 20.016847, 218.0029286))

    assert_allclose(vdist(10, 20, lat2, lon2), (sr, az, a21))


def test_compare_vicenty():
    pyproj = pytest.importorskip('pyproj')
    
    lat2, lon2, a21 = vreckon(10, 20, sr, az)

    p4lon, p4lat, p4a21 = pyproj.Geod(ellps='WGS84').fwd(lon2, lat2, az, sr)
    assert_allclose((p4lon, p4lat, p4a21 % 360.), (lon2, lat2, a21), rtol=0.0025)

    p4az, p4a21, p4sr = pyproj.Geod(ellps='WGS84').inv(20, 10, lon2, lat2)
    assert_allclose((p4az, p4a21 % 360., p4sr), (az, a21, sr))


def test_compare_geodetic():
    pyproj = pytest.importorskip('pyproj')
    
    xyz = pm.geodetic2ecef(*lla0)

    ecef = pyproj.Proj(proj='geocent', ellps='WGS84', datum='WGS84')
    lla = pyproj.Proj(proj='latlong', ellps='WGS84', datum='WGS84')

    assert_allclose(pyproj.transform(lla, ecef, lla0[1], lla0[0], lla0[2]), xyz)
    assert_allclose(pyproj.transform(ecef, lla, *xyz),
                    (lla0[1], lla0[0], lla0[2]))


if __name__ == '__main__':
    pytest.main(['-x', __file__])
