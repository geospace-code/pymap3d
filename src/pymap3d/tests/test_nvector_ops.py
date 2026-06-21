"""Tests for n-vector spherical geometry operations."""

import math

import pymap3d as pm
from pymap3d.nvector import (
    nvector_distance,
    nvector_interpolate,
    nvector_mean,
    nvector_cross_track_distance,
    nvector_intersection,
    geodetic2nvector,
)
import pytest
from pytest import approx


# Mean Earth radius (equal-volume sphere)
R = 6371000.8


def _g2n(lat, lon):
    """Shorthand: geodetic degrees to n-vector."""
    return geodetic2nvector(lat, lon, deg=True)


class TestNvectorDistance:
    def test_same_point(self):
        n = _g2n(45, 90)
        assert nvector_distance(*n, *n) == approx(0, abs=1e-10)

    def test_antipodal(self):
        n1 = _g2n(0, 0)
        n2 = _g2n(0, 180)
        assert nvector_distance(*n1, *n2) == approx(math.pi * R, rel=1e-6)

    def test_pole_to_equator(self):
        n1 = _g2n(90, 0)
        n2 = _g2n(0, 0)
        assert nvector_distance(*n1, *n2) == approx(math.pi / 2 * R, rel=1e-6)

    def test_known_distance(self):
        # London (51.5, -0.1) to Paris (48.9, 2.35) ~338 km (spherical)
        n1 = _g2n(51.5, -0.1)
        n2 = _g2n(48.9, 2.35)
        dist = nvector_distance(*n1, *n2)
        assert dist == approx(337582, rel=0.001)

    def test_custom_radius(self):
        n1 = _g2n(0, 0)
        n2 = _g2n(0, 90)
        dist = nvector_distance(*n1, *n2, radius=1.0)
        assert dist == approx(math.pi / 2, rel=1e-10)


class TestNvectorInterpolate:
    def test_fraction_zero(self):
        n1 = _g2n(0, 0)
        n2 = _g2n(0, 90)
        result = nvector_interpolate(*n1, *n2, 0.0)
        assert result[0] == approx(n1[0], abs=1e-10)
        assert result[1] == approx(n1[1], abs=1e-10)
        assert result[2] == approx(n1[2], abs=1e-10)

    def test_fraction_one(self):
        n1 = _g2n(0, 0)
        n2 = _g2n(0, 90)
        result = nvector_interpolate(*n1, *n2, 1.0)
        assert result[0] == approx(n2[0], abs=1e-10)
        assert result[1] == approx(n2[1], abs=1e-10)
        assert result[2] == approx(n2[2], abs=1e-10)

    def test_midpoint_equator(self):
        # Midpoint between (0,0) and (0,90) should be (0,45)
        n1 = _g2n(0, 0)
        n2 = _g2n(0, 90)
        result = nvector_interpolate(*n1, *n2, 0.5)
        expected = _g2n(0, 45)
        assert result[0] == approx(expected[0], abs=1e-10)
        assert result[1] == approx(expected[1], abs=1e-10)
        assert result[2] == approx(expected[2], abs=1e-10)

    def test_midpoint_meridian(self):
        # Midpoint between (0,0) and (90,0) should be (45,0)
        n1 = _g2n(0, 0)
        n2 = _g2n(90, 0)
        result = nvector_interpolate(*n1, *n2, 0.5)
        expected = _g2n(45, 0)
        assert result[0] == approx(expected[0], abs=1e-10)
        assert result[1] == approx(expected[1], abs=1e-10)
        assert result[2] == approx(expected[2], abs=1e-10)

    def test_same_point(self):
        n = _g2n(30, 60)
        result = nvector_interpolate(*n, *n, 0.5)
        assert result[0] == approx(n[0], abs=1e-10)
        assert result[1] == approx(n[1], abs=1e-10)
        assert result[2] == approx(n[2], abs=1e-10)


class TestNvectorMean:
    def test_single_point(self):
        n = _g2n(45, 90)
        result = nvector_mean([n[0]], [n[1]], [n[2]])
        assert result[0] == approx(n[0], abs=1e-10)
        assert result[1] == approx(n[1], abs=1e-10)
        assert result[2] == approx(n[2], abs=1e-10)

    def test_symmetric_points(self):
        # Points symmetric about the equator at lon=0 should average to (0, 0)
        n1 = _g2n(30, 0)
        n2 = _g2n(-30, 0)
        result = nvector_mean([n1[0], n2[0]], [n1[1], n2[1]], [n1[2], n2[2]])
        expected = _g2n(0, 0)
        assert result[0] == approx(expected[0], abs=1e-10)
        assert result[1] == approx(expected[1], abs=1e-10)
        assert result[2] == approx(expected[2], abs=1e-10)

    def test_four_symmetric_points(self):
        # Four points at (45,0), (45,90), (45,180), (45,270) -> mean at north pole
        pts = [_g2n(45, lon) for lon in [0, 90, 180, 270]]
        n1s = [p[0] for p in pts]
        n2s = [p[1] for p in pts]
        n3s = [p[2] for p in pts]
        result = nvector_mean(n1s, n2s, n3s)
        expected = _g2n(90, 0)
        assert result[2] == approx(expected[2], abs=1e-6)


class TestNvectorCrossTrackDistance:
    def test_point_on_path(self):
        # Point on the equator, path along the equator -> distance = 0
        n1 = _g2n(0, 0)
        n2 = _g2n(0, 90)
        np_ = _g2n(0, 45)
        assert nvector_cross_track_distance(*n1, *n2, *np_) == approx(0, abs=1)

    def test_pole_to_equator_path(self):
        # North pole, path along the equator -> distance = R * pi/2
        n1 = _g2n(0, 0)
        n2 = _g2n(0, 90)
        np_ = _g2n(90, 0)
        dist = nvector_cross_track_distance(*n1, *n2, *np_)
        assert abs(dist) == approx(math.pi / 2 * R, rel=1e-6)

    def test_sign(self):
        # Point north of eastward equatorial path should be positive (left)
        n1 = _g2n(0, 0)
        n2 = _g2n(0, 90)
        np_north = _g2n(10, 45)
        np_south = _g2n(-10, 45)
        d_north = nvector_cross_track_distance(*n1, *n2, *np_north)
        d_south = nvector_cross_track_distance(*n1, *n2, *np_south)
        assert d_north > 0
        assert d_south < 0


class TestNvectorIntersection:
    def test_equator_and_meridian(self):
        # Equator path (0,0)-(0,90) intersects meridian (90,0)-(-90,0)
        # Intersection should be at (0,0) or (0,180)
        n1 = _g2n(0, -10)
        n2 = _g2n(0, 10)
        n3 = _g2n(10, 0)
        n4 = _g2n(-10, 0)
        result = nvector_intersection(*n1, *n2, *n3, *n4)
        # Should be near (0, 0)
        expected = _g2n(0, 0)
        assert result[0] == approx(expected[0], abs=1e-10)
        assert result[1] == approx(expected[1], abs=1e-10)
        assert result[2] == approx(expected[2], abs=1e-10)

    def test_perpendicular_meridians(self):
        # Meridian at lon=0 and meridian at lon=90 intersect at the poles
        n1 = _g2n(10, 0)
        n2 = _g2n(20, 0)
        n3 = _g2n(10, 90)
        n4 = _g2n(20, 90)
        result = nvector_intersection(*n1, *n2, *n3, *n4)
        # Should be near north pole (n3 component ~1)
        assert abs(result[2]) == approx(1.0, abs=1e-10)
