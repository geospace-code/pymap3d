#!/usr/bin/env python

# Copyright (c) 2014-2018 Michael Hirsch, Ph.D.
#
# Redistribution and use in source and binary forms, with or without modification, are permitted provided that the following conditions are met:
# 1. Redistributions of source code must retain the above copyright notice, this list of conditions and the following disclaimer.
# 2. Redistributions in binary form must reproduce the above copyright notice, this list of conditions and the following disclaimer in the documentation and/or other materials provided with the distribution.
# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

from __future__ import division
from numpy import atleast_1d, empty_like, pi, nan
from datetime import datetime
try:
    from astropy.time import Time
    import astropy.units as u
    from astropy.coordinates import Longitude
except ImportError:
    Time = None

#
from .timeconv import str2dt
"""
The "usevallado" datetime to julian runs 4 times faster than astropy.
However, AstroPy is more accurate.

"""


def datetime2sidereal(t, lon_radians, usevallado=True):
    """
    Convert ``datetime`` to sidereal time

    :algorithm: D. Vallado Fundamentals of Astrodynamics and Applications


    t
        Python datetime

    lon
        longitude in RADIANS
    """
    if usevallado:
        jd = datetime2julian(t)
# %% Greenwich Sidereal time RADIANS
        gst = julian2sidereal(jd)
# %% Algorithm 15 p. 188 rotate GST to LOCAL SIDEREAL TIME
        tsr = gst + lon_radians  # radians
    else:  # astropy
        if Time is not None:
            tsr = Time(t).sidereal_time(kind='apparent',
                                    longitude=Longitude(lon_radians, unit=u.radian)).radian
        else:
            raise ImportError('AstroPy required, or use "usevallado=True"')

    return tsr


def datetime2julian(t):
    """
    Python datetime to Julian time

    from D.Vallado Fundamentals of Astrodynamics and Applications p.187
     and J. Meeus Astronomical Algorithms 1991 Eqn. 7.1 pg. 61
    """

    t = str2dt(t)
    t = atleast_1d(t)

    assert isinstance(t[0], datetime)

    jDate = empty_like(t, dtype=float)  # yes we need the dtype!

    for i, d in enumerate(t):
        if d is None:
            jDate[i] = nan
            continue

        if d.month < 3:
            year = d.year - 1
            month = d.month + 12
        else:
            year = d.year
            month = d.month

        A = int(year / 100.0)
        B = 2 - A + int(A / 4.)
        C = ((d.second / 60. + d.minute) / 60. + d.hour) / 24.
        jDate[i] = (int(365.25 * (year + 4716)) +
                    int(30.6001 * (month + 1)) + d.day + B - 1524.5 + C)

    return jDate


def julian2sidereal(juliandate):
    """
    Convert Julian time to sidereal time

    D. Vallado Ed. 4

    input:

    juliandate
        Julian centuries from J2000.0

    """

# %% Vallado Eq. 3-42 p. 184, Seidelmann 3.311-1
    tUT1 = (juliandate - 2451545.0) / 36525.

    gmst_sec = (67310.54841 + (876600 * 3600 + 8640184.812866) *
        tUT1 + 0.093104 * tUT1**2 - 6.2e-6 * tUT1**3) # Eqn. 3-47 p. 188

    # 1/86400 and %(2*pi) implied by units of radians
    return gmst_sec * (2 * pi) / 86400. % (2 * pi)
