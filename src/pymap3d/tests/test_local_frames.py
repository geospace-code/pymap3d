from math import radians

import pymap3d as pm
from pytest import approx

lla0 = (42, -82, 200)
aer0 = (33, 70, 1000)


def _matvec3(matrix, vector):
    return tuple(sum(matrix[i][j] * vector[j] for j in range(3)) for i in range(3))


def _matmul3(a, b):
    return tuple(
        tuple(sum(a[i][k] * b[k][j] for k in range(3)) for j in range(3)) for i in range(3)
    )


def _assert_matrix_close(actual, expected):
    for row_actual, row_expected in zip(actual, expected):
        assert row_actual == approx(row_expected)


def test_sez_conversions():
    enu = pm.aer2enu(*aer0)
    sez = (-enu[1], enu[0], enu[2])
    xyz = pm.aer2ecef(*aer0, *lla0)

    assert pm.aer2sez(*aer0) == approx(sez)
    assert pm.sez2aer(*sez) == approx(aer0)
    assert pm.ecef2sez(*xyz, *lla0) == approx(sez)
    assert pm.sez2ecef(*sez, *lla0) == approx(xyz)

    lla1 = pm.aer2geodetic(*aer0, *lla0)
    assert pm.geodetic2sez(*lla1, *lla0) == approx(sez)
    assert pm.sez2geodetic(*sez, *lla0) == approx(lla1)


def test_nwu_conversions():
    enu = pm.aer2enu(*aer0)
    nwu = (enu[1], -enu[0], enu[2])
    xyz = pm.aer2ecef(*aer0, *lla0)

    assert pm.aer2nwu(*aer0) == approx(nwu)
    assert pm.nwu2aer(*nwu) == approx(aer0)
    assert pm.ecef2nwu(*xyz, *lla0) == approx(nwu)
    assert pm.nwu2ecef(*nwu, *lla0) == approx(xyz)

    lla1 = pm.aer2geodetic(*aer0, *lla0)
    assert pm.geodetic2nwu(*lla1, *lla0) == approx(nwu)
    assert pm.nwu2geodetic(*nwu, *lla0) == approx(lla1)


def test_vector_frame_conversions():
    ecef_vector = (5, 3, 2)
    enu_vector = pm.ecef2enuv(*ecef_vector, *lla0[:2])
    sez_vector = (-enu_vector[1], enu_vector[0], enu_vector[2])
    nwu_vector = (enu_vector[1], -enu_vector[0], enu_vector[2])

    assert pm.ecef2sezv(*ecef_vector, *lla0[:2]) == approx(sez_vector)
    assert pm.sez2ecefv(*sez_vector, *lla0[:2]) == approx(ecef_vector)
    assert pm.ecef2nwuv(*ecef_vector, *lla0[:2]) == approx(nwu_vector)
    assert pm.nwu2ecefv(*nwu_vector, *lla0[:2]) == approx(ecef_vector)


def test_frame_matrices_match_vector_transforms():
    ecef_vector = (5, 3, 2)
    lat0, lon0 = lla0[:2]

    frame_cases = [
        (pm.ecef2enu_matrix(lat0, lon0), pm.enu2ecef_matrix(lat0, lon0), pm.ecef2enuv(*ecef_vector, lat0, lon0)),
        (pm.ecef2ned_matrix(lat0, lon0), pm.ned2ecef_matrix(lat0, lon0), pm.ecef2nedv(*ecef_vector, lat0, lon0)),
        (pm.ecef2sez_matrix(lat0, lon0), pm.sez2ecef_matrix(lat0, lon0), pm.ecef2sezv(*ecef_vector, lat0, lon0)),
        (pm.ecef2nwu_matrix(lat0, lon0), pm.nwu2ecef_matrix(lat0, lon0), pm.ecef2nwuv(*ecef_vector, lat0, lon0)),
    ]

    identity = ((1.0, 0.0, 0.0), (0.0, 1.0, 0.0), (0.0, 0.0, 1.0))

    for forward, inverse, expected in frame_cases:
        assert _matvec3(forward, ecef_vector) == approx(expected)
        assert _matvec3(inverse, expected) == approx(ecef_vector)
        _assert_matrix_close(_matmul3(forward, inverse), identity)


def test_frame_matrices_rad_inputs():
    lat0, lon0 = lla0[:2]
    rlat0 = radians(lat0)
    rlon0 = radians(lon0)

    _assert_matrix_close(pm.ecef2enu_matrix(lat0, lon0), pm.ecef2enu_matrix(rlat0, rlon0, deg=False))
    _assert_matrix_close(pm.ecef2ned_matrix(lat0, lon0), pm.ecef2ned_matrix(rlat0, rlon0, deg=False))
    _assert_matrix_close(pm.ecef2sez_matrix(lat0, lon0), pm.ecef2sez_matrix(rlat0, rlon0, deg=False))
    _assert_matrix_close(pm.ecef2nwu_matrix(lat0, lon0), pm.ecef2nwu_matrix(rlat0, rlon0, deg=False))
