#!/usr/bin/env python
from numpy import cos, arcsin, sqrt, radians, degrees
from astropy.coordinates.angle_utilities import angular_separation
"""
Michael Hirsch


inputs:
r0,d0: for first point, rightAscension,Declination [degrees]  or (azimuth,elevation)
r1,d1: for second point, rightAscension,Declination [degrees] or (azimuth,elevation)

(or, azimuth/elevation respectively)

Note: adding decimal points to the constants made 0 difference in %timeit execution time

The Meeus algorithm is about 9.5% faster than Astropy/Vicenty on my PC,
and gives virtually identical result
within double precision arithmetic limitations
"""


def angledist_meeus(r0, d0, r1, d1):
    """
    Meeus
    from "Astronomical Algorithms" by Jean Meeus Ch. 16 p. 111 (16.5)
    gives angular distance in degrees between two rightAscension,Declination
    points in the sky.  Neglecting atmospheric effects, of course.

    Advantage of Meeus haversine method is stability all the way to exactly 0 deg.

    assumes degrees input, degrees output

    either the arrays must be the same size, or one of them must be a scalar
    """

    r0 = radians(r0)
    r1 = radians(r1)
    d0 = radians(d0)
    d1 = radians(d1)
    dist_rad = 2 * arcsin(
        sqrt(
            haversine(d1 - d0) +
            cos(d0) * cos(d1) * haversine(r1 - r0)))

    return degrees(dist_rad)


def angledist(lon1, lat1, lon2, lat2):
    """
    For reference, this is from astropy astropy/coordinates/angle_utilities.py
    Angular separation between two points on a sphere.
    """
    return degrees(angular_separation(radians(lon1),
                                      radians(lat1), radians(lon2),
                                      radians(lat2)))


def haversine(theta):
    """
    http://en.wikipedia.org/wiki/Haversine
    Meeus p. 111 """
    return (1 - cos(theta)) / 2.
