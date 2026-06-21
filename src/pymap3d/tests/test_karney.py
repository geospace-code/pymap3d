"""Tests for Karney's geodesic algorithms."""

import math

import pymap3d.karney as karney
import pytest
from pytest import approx


class TestGeodesicInverse:
    def test_coincident_points(self):
        dist, azi1, azi2 = karney.geodesic_inverse(0, 0, 0, 0)
        assert dist == approx(0, abs=1e-6)

    def test_equatorial(self):
        """Quarter circumference along equator."""
        dist, azi1, azi2 = karney.geodesic_inverse(0, 0, 0, 90)
        assert dist == approx(1.001875e7, rel=0.001)
        assert azi1 == approx(90, abs=0.01)

    def test_meridional(self):
        """Pole to equator along meridian."""
        dist, azi1, azi2 = karney.geodesic_inverse(90, 0, 0, 0)
        assert dist == approx(1.00019657e7, rel=0.001)
        assert azi1 == approx(180, abs=0.01)

    def test_pole_to_pole(self):
        """North pole to south pole."""
        dist, azi1, azi2 = karney.geodesic_inverse(90, 0, -90, 0)
        assert dist == approx(2.000393145e7, rel=0.001)

    def test_near_antipodal(self):
        """Near-antipodal case where Vincenty struggles."""
        dist, azi1, azi2 = karney.geodesic_inverse(0, 0, 0.5, 179.5)
        assert dist == approx(1.9936e7, rel=0.01)

    def test_short_distance(self):
        """Short distance (~3km)."""
        dist, azi1, azi2 = karney.geodesic_inverse(10, 20, 10.02137267, 20.0168471)
        assert dist == approx(3e3, rel=0.01)

    def test_radians(self):
        dist, azi1, azi2 = karney.geodesic_inverse(0, 0, 0, math.pi / 2, deg=False)
        assert dist == approx(1.001875e7, rel=0.001)

    def test_agreement_with_vincenty(self):
        """For non-pathological cases, should agree with Vincenty."""
        pytest.importorskip("pymap3d.vincenty")
        import pymap3d.vincenty as vincenty

        lat1, lon1 = 40.0, -74.0
        lat2, lon2 = 48.8, 2.35

        dist_k, azi1_k, _ = karney.geodesic_inverse(lat1, lon1, lat2, lon2)
        dist_v, azi1_v = vincenty.vdist(lat1, lon1, lat2, lon2)

        assert dist_k == approx(dist_v, rel=1e-5)
        assert azi1_k == approx(azi1_v, rel=1e-4)


class TestGeodesicDirect:
    def test_zero_distance(self):
        lat2, lon2, azi2 = karney.geodesic_direct(0, 0, 90, 0)
        assert lat2 == approx(0, abs=1e-10)
        assert lon2 == approx(0, abs=1e-10)

    def test_equatorial_east(self):
        """Travel east along equator."""
        lat2, lon2, azi2 = karney.geodesic_direct(0, 0, 90, 1.001875e7)
        assert lat2 == approx(0, abs=0.01)
        assert lon2 == approx(90, abs=0.1)

    def test_north_from_equator(self):
        """Travel north from equator."""
        lat2, lon2, azi2 = karney.geodesic_direct(0, 0, 0, 1.00019657e7)
        assert lat2 == approx(90, abs=0.1)

    def test_roundtrip(self):
        """Direct then inverse should recover distance and azimuth."""
        lat1, lon1, azi1 = 40.0, -74.0, 45.0
        d = 1e6

        lat2, lon2, azi2 = karney.geodesic_direct(lat1, lon1, azi1, d)
        dist_back, azi1_back, _ = karney.geodesic_inverse(lat1, lon1, lat2, lon2)

        assert dist_back == approx(d, rel=1e-6)
        assert azi1_back == approx(azi1, abs=0.01)


class TestGeodesicLine:
    def test_endpoints(self):
        """First and last points should match the given endpoints."""
        lats, lons = karney.geodesic_line(0, 0, 0, 90, npts=10)
        assert lats[0] == approx(0, abs=1e-6)
        assert lons[0] == approx(0, abs=1e-6)
        assert lats[-1] == approx(0, abs=0.1)
        assert lons[-1] == approx(90, abs=0.1)

    def test_npts(self):
        lats, lons = karney.geodesic_line(0, 0, 45, 45, npts=50)
        assert len(lats) == 50
        assert len(lons) == 50


class TestGeodesicArea:
    def test_triangle(self):
        """Small equilateral-ish triangle."""
        lats = [0, 1, 0.5]
        lons = [0, 0, 0.866]
        area, perim = karney.geodesic_area(lats, lons)
        # Should be positive (CCW) and roughly match planar area
        assert area > 0
        assert perim > 0

    def test_too_few_points(self):
        area, perim = karney.geodesic_area([0, 1], [0, 1])
        assert area == 0.0
        assert perim == 0.0
