#!/usr/bin/env python3
from __future__ import division
from numpy import cos,arcsin,sqrt,radians,degrees,nan
"""
Michael Hirsch
from "Astronomical Algorithms" by Jean Meeus Ch. 16 p. 111 (16.5)
gives angular distance in degrees between two rightAscension,Declination 
points in the sky

inputs:
r0,d0: for first point, rightAscension,Declination [degrees]
r1,d1: for second point, rightAscension,Declination [degrees]

Note: adding decimal points to the constants made 0 difference in 
%timeit execution time
GPLv3+ license
"""

def angledist(r0,d0,r1,d1):
    #assumes degrees input, degrees output
    r0 = radians(r0); r1 = radians(r1)
    d0 = radians(d0); d1 = radians(d1)
    dist_rad = 2 * arcsin( sqrt( haversine(d1-d0) + cos(d0)*cos(d1) * 
haversine(r1-r0) ) )

    return degrees(dist_rad)

def haversine(theta):
    """ http://en.wikipedia.org/wiki/Haversine and Meeus p. 111 """
    return (1-cos(theta)) / 2

if __name__ == '__main__': #pragma: no cover
    from argparse import ArgumentParser
    p = ArgumentParser(description="computes angular distance between 
two points in sky")
    p.add_argument('r0',help='right ascension of first point 
[degrees]',type=float,nargs='?',default=nan)
    p.add_argument('d0',help='declination of first point 
[degrees]',type=float,nargs='?',default=nan)
    p.add_argument('r1',help='right ascension of second point 
[degrees]',type=float,nargs='?',default=nan)
    p.add_argument('d1',help='declination of second point 
[degrees]',type=float,nargs='?',default=nan)
    a = p.parse_args()
    
    dist_deg = angledist(a.r0,a.d0,a.r1,a.d1)
    print(dist_deg)


