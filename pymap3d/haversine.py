#!/usr/bin/env python3
from __future__ import division
from numpy import cos,sin,arctan2,arcsin,sqrt,radians,degrees
"""
Michael Hirsch
from "Astronomical Algorithms" by Jean Meeus Ch. 16 p. 111 (16.5)
gives angular distance in degrees between two rightAscension,Declination
points in the sky.  Neglecting atmospheric effects, of course.

inputs:
r0,d0: for first point, rightAscension,Declination [degrees]  or (azimuth,elevation)
r1,d1: for second point, rightAscension,Declination [degrees] or (azimuth,elevation)

(or, azimuth/elevation respectively)

Note: adding decimal points to the constants made 0 difference in %timeit execution time
GPLv3+ license

The Meeus algorithm is about 9.5% faster than Astropy/Vicenty on my PC, and gives virtually identical result
within double precision arithmetic limitations
"""

def angledist(r0,d0,r1,d1):
    """ Meeus
    assumes degrees input, degrees output
    """

    r0 = radians(r0); r1 = radians(r1)
    d0 = radians(d0); d1 = radians(d1)
    dist_rad = 2*arcsin(
                 sqrt(
                 haversine(d1-d0) +
                   cos(d0)*cos(d1) * haversine(r1-r0) ) )

    return degrees(dist_rad)

def angular_separation(lon1, lat1, lon2, lat2):
    """
    For reference, this is from astropy astropy/coordinates/angle_utilities.py
    Angular separation between two points on a sphere.
    Parameters
    ----------
    lon1, lat1, lon2, lat2 : `Angle`, `~astropy.units.Quantity` or float
        Longitude and latitude of the two points. Quantities should be in
        angular units; floats in radians.
    Returns
    -------
    angular separation : `~astropy.units.Quantity` or float
        Type depends on input; `Quantity` in angular units, or float in
        radians.
    Notes
    -----
    The angular separation is calculated using the Vincenty formula [1]_,
    which is slightly more complex and computationally expensive than
    some alternatives, but is stable at at all distances, including the
    poles and antipodes.
    .. [1] http://en.wikipedia.org/wiki/Great-circle_distance
    """
    sdlon = sin(lon2 - lon1)
    cdlon = cos(lon2 - lon1)
    slat1 = sin(lat1)
    slat2 = sin(lat2)
    clat1 = cos(lat1)
    clat2 = cos(lat2)

    num1 = clat2 * sdlon
    num2 = clat1 * slat2 - slat1 * clat2 * cdlon
    denominator = slat1 * slat2 + clat1 * clat2 * cdlon

    return degrees(arctan2(sqrt(num1 ** 2 + num2 ** 2), denominator))

def haversine(theta):
    """
    http://en.wikipedia.org/wiki/Haversine
    Meeus p. 111 """
    return (1-cos(theta)) / 2

if __name__ == '__main__': #pragma: no cover
    from argparse import ArgumentParser
    p = ArgumentParser(description="computes angular distance between two points in sky")
    p.add_argument('r0',help='right ascension of first point [degrees]',type=float)
    p.add_argument('d0',help='declination of first point [degrees]',type=float)
    p.add_argument('r1',help='right ascension of second point [degrees]',type=float)
    p.add_argument('d1',help='declination of second point [degrees]', type=float)
    a = p.parse_args()

    dist_deg = angledist(a.r0,a.d0,a.r1,a.d1)

    print(dist_deg)
