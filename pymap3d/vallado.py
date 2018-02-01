#!/usr/bin/env python

# Copyright (c) 2014-2018 Michael Hirsch, Ph.D.
#
# Redistribution and use in source and binary forms, with or without modification, are permitted provided that the following conditions are met:
# 1. Redistributions of source code must retain the above copyright notice, this list of conditions and the following disclaimer.
# 2. Redistributions in binary form must reproduce the above copyright notice, this list of conditions and the following disclaimer in the documentation and/or other materials provided with the distribution.
# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

"""
converts right ascension, declination to azimuth, elevation and vice versa.
Normally do this via AstroPy.
These functions are fallbacks for those who don't wish to use AstroPy (perhaps Python 2.7 users).

Michael Hirsch implementation of algorithms from D. Vallado
"""
from __future__ import division
try:
    import numpy
    from numpy import sin, cos, degrees, radians, arcsin, arctan2, atleast_1d
except ImportError:
    numpy = None
    from math import sin, cos, degrees, radians
    from math import atan2 as arctan2
    from math import asin as arcsin
#
from .datetime2hourangle import datetime2sidereal


def vazel2radec(az_deg, el_deg, lat_deg, lon_deg, t):
    """
    convert azimuth, elevation to right ascension, declination

    Inputs

    az_deg
        Numpy ndarray of azimuth to point [degrees]

    el_deg
        Numpy ndarray of elevation to point [degrees]

    lat_deg
        scalar observer WGS84 latitude [degrees]

    lon_deg
        scalar observer WGS84 longitude [degrees]

    t
        time of observation

    Outputs

    ra_deg
        Numpy ndarray of right ascension values [degrees]

    dec_deg
        Numpy ndarray of declination values [degrees]

    from D.Vallado Fundamentals of Astrodynamics and Applications
    p.258-259
    """
    if numpy is not None:
        az_deg = atleast_1d(az_deg)
        el_deg = atleast_1d(el_deg)
        lat_deg = atleast_1d(lat_deg)
        lon_deg = atleast_1d(lon_deg)

        assert az_deg.shape == el_deg.shape, 'az and el must be same shape ndarray'
        assert lat_deg.size == 1 and lon_deg.size == 1, 'need one observer and one or more  (az,el).'


    az = radians(az_deg)
    el = radians(el_deg)
    lat = radians(lat_deg)
    lon = radians(lon_deg)
# %% Vallado "algorithm 28" p 268
    dec = arcsin(sin(el) * sin(lat) + cos(el) * cos(lat) * cos(az))

    lha = arctan2(-(sin(az) * cos(el)) / cos(dec),
                  (sin(el) - sin(lat) * sin(dec)) / (cos(dec) * cos(lat)))

    lst = datetime2sidereal(t, lon)  # lon, ra in RADIANS

    """ by definition right ascension \in [0,360) degrees """
    return degrees(lst - lha) % 360, degrees(dec)


def vradec2azel(ra_deg,dec_deg,lat_deg,lon_deg,t):
    """
    convert right ascension, declination to azimuth, elevation

    Inputs

    ra_deg
        Numpy ndarray of right ascension values [degrees]

    dec_deg
        Numpy ndarray of declination values [degrees]

    lat_deg
        scalar observer WGS84 latitude [degrees]

    lon_deg
        scalar observer WGS84 longitude [degrees]

    t
        time of observation

    Outputs

    az_deg
        Numpy ndarray of azimuth to point [degrees]

    el_deg
        Numpy ndarray of elevation to point [degrees]


    from D. Vallado "Fundamentals of Astrodynamics and Applications "
       4th Edition Ch. 4.4 pg. 266-268
    """
    ra =  radians(ra_deg)
    dec = radians(dec_deg)
    lat = radians(lat_deg)
    lon = radians(lon_deg)

    lst = datetime2sidereal(t, lon) #RADIANS
# %% Eq. 4-11 p. 267 LOCAL HOUR ANGLE
    lha = lst - ra
# %% #Eq. 4-12 p. 267
    el = arcsin(sin(lat) * sin(dec) + cos(lat) * cos(dec) * cos(lha) )
# %% combine Eq. 4-13 and 4-14 p. 268
    az = arctan2( -sin(lha) * cos(dec) / cos(el),
                   ( sin(dec) - sin(el) * sin(lat) ) / (cos(el)*cos(lat))  )


    return degrees(az) % 360.0, degrees(el)
