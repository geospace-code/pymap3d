from __future__ import annotations

from typing import Any

import pytest
from pytest import approx

try:
    from numpy import asarray
    from numpy.typing import NDArray
except ImportError:
    pass


import pymap3d as pm

ELL = pm.Ellipsoid.from_name("wgs84")
A = ELL.semimajor_axis
B = ELL.semiminor_axis

llrlla: list[tuple[tuple[float, float, float], tuple[float, float, float]]] = [
    ((0, 0, A - 1), (0, 0, -1)),
    ((0, 90, A - 1), (0, 90, -1)),
    ((0, -90, A + 1), (0, -90, 1)),
    ((44.807576814237606, 270, 6367490.543857), (45, 270, 1)),
    ((90, 0, B + 1), (90, 0, 1)),
    ((90, 15, B - 1), (90, 15, -1)),
    ((-90, 0, B + 1), (-90, 0, 1)),
]
llallr: list[tuple[tuple[float, float, float], tuple[float, float, float]]] = [
    ((0, 0, -1), (0, 0, A - 1)),
    ((0, 90, -1), (0, 90, A - 1)),
    ((0, -90, 1), (0, -90, A + 1)),
    ((45, 270, 1), (44.807576814237606, 270, 6367490.543857)),
    ((90, 0, 1), (90, 0, B + 1)),
    ((90, 15, -1), (90, 15, B - 1)),
    ((-90, 0, 1), (-90, 0, B + 1)),
]
llallr_list = [([[i] for i in lla], llr) for lla, llr in llallr]
llrlla_list = [([[i] for i in llr], lla) for llr, lla in llrlla]
try:
    llallr_array = [([asarray(i) for i in lla], llr) for lla, llr in llallr]
    llrlla_array = [([asarray(i) for i in llr], lla) for llr, lla in llrlla]
except NameError:
    llallr_array = []
    llrlla_array = []


atol_dist = 1e-6  # 1 micrometer


@pytest.mark.parametrize("lla, llr", llallr)
def test_geodetic2spherical(
    lla: tuple[float, float, float], llr: tuple[float, float, float]
) -> None:
    coords = pm.geodetic2spherical(*lla)
    assert coords[:2] == approx(llr[:2])
    assert coords[2] == approx(llr[2], abs=atol_dist)


@pytest.mark.parametrize("llr, lla", llrlla)
def test_spherical2geodetic(
    llr: tuple[float, float, float], lla: tuple[float, float, float]
) -> None:
    coords = pm.spherical2geodetic(*llr)
    assert coords[:2] == approx(lla[:2])
    assert coords[2] == approx(lla[2], abs=atol_dist)


@pytest.mark.parametrize("lla, llr", llallr_list)
def test_geodetic2spherical_list(lla: list[list[float]], llr: tuple[float, float, float]) -> None:
    pytest.importorskip("numpy")
    coords = pm.geodetic2spherical(*lla)  # type: ignore[call-overload]
    assert coords[:2] == approx(llr[:2])
    assert coords[2] == approx(llr[2], abs=atol_dist)


@pytest.mark.parametrize("llr, lla", llrlla_list)
def test_spherical2geodetic_list(llr: list[list[float]], lla: tuple[float, float, float]) -> None:
    pytest.importorskip("numpy")
    coords = pm.spherical2geodetic(*llr)  # type: ignore[call-overload]
    assert coords[:2] == approx(lla[:2])
    assert coords[2] == approx(lla[2], abs=atol_dist)


@pytest.mark.parametrize("lla, llr", llallr_array)
def test_geodetic2spherical_array(lla: list[NDArray[Any]], llr: tuple[float, float, float]) -> None:
    pytest.importorskip("numpy")
    coords = pm.geodetic2spherical(*lla)  # type: ignore[call-overload]
    assert coords[:2] == approx(llr[:2])
    assert coords[2] == approx(llr[2], abs=atol_dist)


@pytest.mark.parametrize("llr, lla", llrlla_array)
def test_spherical2geodetic_array(llr: list[NDArray[Any]], lla: tuple[float, float, float]) -> None:
    pytest.importorskip("numpy")
    coords = pm.spherical2geodetic(*llr)  # type: ignore[call-overload]
    assert coords[:2] == approx(lla[:2])
    assert coords[2] == approx(lla[2], abs=atol_dist)
