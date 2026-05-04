from math import radians

import pymap3d as pm
import pytest
from pytest import approx

lla0 = (42, -82, 200)
ypr0 = (30, 12, -8)


def _matvec3(matrix, vector):
    return tuple(sum(matrix[i][j] * vector[j] for j in range(3)) for i in range(3))


def _matmul3(a, b):
    return tuple(
        tuple(sum(a[i][k] * b[k][j] for k in range(3)) for j in range(3)) for i in range(3)
    )


def _assert_matrix_close(actual, expected):
    for row_actual, row_expected in zip(actual, expected):
        assert row_actual == approx(row_expected)


def test_body_matrix_identity():
    identity = ((1.0, 0.0, 0.0), (0.0, 1.0, 0.0), (0.0, 0.0, 1.0))

    _assert_matrix_close(pm.body_matrix(0, 0, 0), identity)
    _assert_matrix_close(pm.ned2body_matrix(0, 0, 0), identity)
    _assert_matrix_close(pm.body2ned_matrix(0, 0, 0), identity)


def test_body_rotation_conventions():
    assert pm.ned2body(1, 0, 0, 0, 0, 0) == approx((1, 0, 0))
    assert pm.ned2body(0, 1, 0, 90, 0, 0) == approx((1, 0, 0), abs=1e-12)
    assert pm.ned2body(0, 0, -1, 0, 90, 0) == approx((1, 0, 0), abs=1e-12)
    assert pm.ned2body(0, 0, 1, 0, 0, 90) == approx((0, 1, 0), abs=1e-12)


def test_body_matrix_inverse_and_radians():
    yaw, pitch, roll = ypr0
    matrix = pm.ned2body_matrix(yaw, pitch, roll)
    inverse = pm.body2ned_matrix(yaw, pitch, roll)
    identity = ((1.0, 0.0, 0.0), (0.0, 1.0, 0.0), (0.0, 0.0, 1.0))

    _assert_matrix_close(_matmul3(matrix, inverse), identity)
    _assert_matrix_close(
        matrix,
        pm.ned2body_matrix(radians(yaw), radians(pitch), radians(roll), deg=False),
    )


def test_ned_enu_body_consistency():
    enu = (5, 3, 2)
    ned = (enu[1], enu[0], -enu[2])
    sez = (-enu[1], enu[0], enu[2])
    nwu = (enu[1], -enu[0], enu[2])
    yaw, pitch, roll = ypr0

    assert pm.enu2body(*enu, yaw, pitch, roll) == approx(pm.ned2body(*ned, yaw, pitch, roll))
    assert pm.sez2body(*sez, yaw, pitch, roll) == approx(pm.enu2body(*enu, yaw, pitch, roll))
    assert pm.nwu2body(*nwu, yaw, pitch, roll) == approx(pm.enu2body(*enu, yaw, pitch, roll))
    assert pm.body2enu(*pm.enu2body(*enu, yaw, pitch, roll), yaw, pitch, roll) == approx(enu)
    assert pm.body2ned(*pm.ned2body(*ned, yaw, pitch, roll), yaw, pitch, roll) == approx(ned)
    assert pm.body2sez(*pm.sez2body(*sez, yaw, pitch, roll), yaw, pitch, roll) == approx(sez)
    assert pm.body2nwu(*pm.nwu2body(*nwu, yaw, pitch, roll), yaw, pitch, roll) == approx(nwu)


def test_body_matrices_match_vector_transforms():
    local_ned = (4, -2, 1)
    local_enu = (local_ned[1], local_ned[0], -local_ned[2])
    local_sez = (-local_enu[1], local_enu[0], local_enu[2])
    local_nwu = (local_enu[1], -local_enu[0], local_enu[2])
    yaw, pitch, roll = ypr0

    assert _matvec3(pm.ned2body_matrix(yaw, pitch, roll), local_ned) == approx(
        pm.ned2body(*local_ned, yaw, pitch, roll)
    )
    assert _matvec3(pm.enu2body_matrix(yaw, pitch, roll), local_enu) == approx(
        pm.enu2body(*local_enu, yaw, pitch, roll)
    )
    assert _matvec3(pm.sez2body_matrix(yaw, pitch, roll), local_sez) == approx(
        pm.sez2body(*local_sez, yaw, pitch, roll)
    )
    assert _matvec3(pm.nwu2body_matrix(yaw, pitch, roll), local_nwu) == approx(
        pm.nwu2body(*local_nwu, yaw, pitch, roll)
    )


@pytest.mark.parametrize("frame", ["ned", "enu", "sez", "nwu"])
def test_ecef_body_roundtrip(frame):
    yaw, pitch, roll = ypr0
    target_ecef = pm.aer2ecef(33, 40, 1200, *lla0)

    body = pm.ecef2body(*target_ecef, *lla0, yaw, pitch, roll, frame=frame)
    assert pm.body2ecef(*body, *lla0, yaw, pitch, roll, frame=frame) == approx(target_ecef)


@pytest.mark.parametrize("frame", ["ned", "enu", "sez", "nwu"])
def test_ecef_body_vector_roundtrip(frame):
    yaw, pitch, roll = ypr0
    ecef_vector = (5, 3, 2)

    body_vector = pm.ecef2bodyv(*ecef_vector, *lla0[:2], yaw, pitch, roll, frame=frame)
    assert pm.body2ecefv(*body_vector, *lla0[:2], yaw, pitch, roll, frame=frame) == approx(
        ecef_vector
    )


@pytest.mark.parametrize("frame", ["ned", "enu", "sez", "nwu"])
def test_geodetic_body_roundtrip(frame):
    yaw, pitch, roll = ypr0
    target_lla = pm.aer2geodetic(33, 40, 1200, *lla0)

    body = pm.geodetic2body(*target_lla, *lla0, yaw, pitch, roll, frame=frame)
    assert pm.body2geodetic(*body, *lla0, yaw, pitch, roll, frame=frame) == approx(target_lla)


def test_body_invalid_options():
    with pytest.raises(ValueError):
        pm.body_matrix(0, 0, 0, sequence="xyz")

    with pytest.raises(ValueError):
        pm.ecef2body(1, 2, 3, *lla0, 0, 0, 0, frame="xyz")
