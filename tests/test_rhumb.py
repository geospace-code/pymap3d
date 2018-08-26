#!/usr/bin/env python
import pytest
from pytest import approx
import pymap3d as pm


def test_meridianarc():
    """
    meridianarc(0, deg2rad(40), wgs84Ellipsoid)
    """

    mdist = pm.meridian_dist(40)

    assert mdist == approx(4429529.03035058)


def test_rhumb():
    """
    distance('rh', 40, -80, 65, -148, wgs84Ellipsoid)
    azimuth('rh', 40, -80, 65, -148, wgs84Ellipsoid)
    """
    rhdist, rhaz = pm.loxodrome_inverse(40, -80, 65, -148)

    assert rhdist == approx(5248666.20853187)
    assert rhaz == approx(302.0056736)


if __name__ == '__main__':
    pytest.main(['-xrsv', __file__])
