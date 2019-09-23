#!/usr/bin/env python
import pytest
from pytest import approx
from math import radians

import pymap3d as pm


@pytest.mark.parametrize("geodetic_lat,geocentric_lat", [(0, 0), (90, 90), (-90, -90), (45, 44.80757678), (-45, -44.80757678)])
def test_geodetic_geocentric(geodetic_lat, geocentric_lat):

    assert pm.geodetic2geocentric(geodetic_lat) == approx(geocentric_lat)
    assert pm.geodetic2geocentric(radians(geodetic_lat), deg=False) == approx(radians(geocentric_lat))

    assert pm.geocentric2geodetic(geocentric_lat) == approx(geodetic_lat)
    assert pm.geocentric2geodetic(radians(geocentric_lat), deg=False) == approx(radians(geodetic_lat))


def test_numpy_geodetic_geocentric():
    pytest.importorskip("numpy")
    assert pm.geodetic2geocentric([45, 0]) == approx([44.80757678, 0])
    assert pm.geocentric2geodetic([44.80757678, 0]) == approx([45, 0])


@pytest.mark.parametrize("lat", [91, -91])
def test_badvals(lat):
    with pytest.raises(ValueError):
        pm.geodetic2geocentric(lat)
        pm.geocentric2geodetic(lat)


if __name__ == "__main__":
    pytest.main(["-v", __file__])
