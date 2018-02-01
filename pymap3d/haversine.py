#!/usr/bin/env python
# Copyright (c) 2014-2018 Michael Hirsch, Ph.D.
#
# Redistribution and use in source and binary forms, with or without modification, are permitted provided that the following conditions are met:
# 1. Redistributions of source code must retain the above copyright notice, this list of conditions and the following disclaimer.
# 2. Redistributions in binary form must reproduce the above copyright notice, this list of conditions and the following disclaimer in the documentation and/or other materials provided with the distribution.
# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

try:
    from numpy import cos, arcsin, sqrt, radians, degrees
except ImportError:
    from math import cos,sqrt,radians,degrees
    from math import asin as arcsin
try:
    from astropy.coordinates.angle_utilities import angular_separation
except ImportError:
    angular_separation = None
"""
Michael Hirsch

Note: adding decimal points to the constants made 0 difference in %timeit execution time

The Meeus algorithm is about 9.5% faster than Astropy/Vicenty on my PC,
and gives virtually identical result
within double precision arithmetic limitations
"""


def angledist_meeus(r0, d0, r1, d1):
    """
    inputs:

    r0,d0
        for first point, rightAscension,Declination [degrees]  or (azimuth,elevation)

    r1,d1
        for second point, rightAscension,Declination [degrees] or (azimuth,elevation)


    (or, azimuth/elevation respectively)


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
    inputs:

    r0,d0
        for first point, rightAscension,Declination [degrees]  or (azimuth,elevation)

    r1,d1
        for second point, rightAscension,Declination [degrees] or (azimuth,elevation)


    (or, azimuth/elevation respectively)

    For reference, this is from astropy astropy/coordinates/angle_utilities.py
    Angular separation between two points on a sphere.
    """
    if angular_separation is None:
        raise ImportError('angledist requires AstroPy. Try pure Python angledis_meeus')

    return degrees(angular_separation(radians(lon1),
                                      radians(lat1), radians(lon2),
                                      radians(lat2)))


def haversine(theta):
    """
    Compute haversine of angle theta (radians)

    http://en.wikipedia.org/wiki/Haversine
    Meeus p. 111
    """
    return (1 - cos(theta)) / 2.
