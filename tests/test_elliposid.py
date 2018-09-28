#!/usr/bin/env python
import pytest
from pytest import approx
import pymap3d as pm

xyz0 = (660e3, -4700e3, 4247e3)


def test_ellipsoid():

    assert pm.ecef2geodetic(*xyz0, ell=pm.Ellipsoid('wgs84')) == approx([42.014670535, -82.0064785, 276.9136916])
    assert pm.ecef2geodetic(*xyz0, ell=pm.Ellipsoid('grs80')) == approx([42.014670536, -82.0064785, 276.9137385])
    assert pm.ecef2geodetic(*xyz0, ell=pm.Ellipsoid('clrk66')) == approx([42.01680003, -82.0064785, 313.9026793])
    assert pm.ecef2geodetic(*xyz0, ell=pm.Ellipsoid('mars')) == approx([42.009428417, -82.006479, 2.981246e6])
    assert pm.ecef2geodetic(*xyz0, ell=pm.Ellipsoid('venus')) == approx([41.8233663, -82.0064785, 3.17878159e5])
    assert pm.ecef2geodetic(*xyz0, ell=pm.Ellipsoid('moon')) == approx([41.8233663, -82.0064785, 4.630878e6])


if __name__ == '__main__':
    pytest.main([__file__])
