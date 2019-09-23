#!/usr/bin/env python
# Copyright (c) 2014-2018 Michael Hirsch, Ph.D.
"""
PyMap3D provides coordinate transforms and geodesy functions with a similar API
to the Matlab Mapping Toolbox, but was of course independently derived.

For all functions, the default units are:

distance : float
    METERS
angles : float
    DEGREES
time : datetime.datetime
    UTC time of observation

These functions may be used with any planetary body, provided the appropriate
reference ellipsoid is defined. The default ellipsoid is WGS-84

deg : bool = True means degrees. False = radians.

Most functions accept NumPy arrays of any shape, as well as compatible data types
including AstroPy, Pandas and Xarray that have Numpy-like data properties.
For clarity, we omit all these types in the docs, and just specify the scalar type.

Other languages
---------------

Companion packages exist for:

* Matlab / GNU Octave: [Matmap3D](https://github.com/scivision/matmap3d)
* Fortran: [Maptran3D](https://github.com/scivision/maptran3d)
"""

from .aer import ecef2aer, aer2ecef, geodetic2aer, aer2geodetic, eci2aer, aer2eci
from .ellipsoid import Ellipsoid
from .enu import enu2geodetic, geodetic2enu, aer2enu, enu2aer
from .ned import ned2ecef, ned2geodetic, geodetic2ned, ecef2nedv, ned2aer, aer2ned, ecef2ned
from .ecef import geodetic2ecef, ecef2geodetic, eci2geodetic, ecef2enuv, enu2ecef, ecef2enu, enu2uvw, uvw2enu
from .sidereal import datetime2sidereal
from .latitude import geodetic2geocentric, geocentric2geodetic
from .lox import isometric, meridian_dist, loxodrome_inverse
from .los import lookAtSpheroid
from .timeconv import str2dt

try:
    from .azelradec import radec2azel, azel2radec
    from .eci import eci2ecef, ecef2eci
except ImportError:
    from .vallado import radec2azel, azel2radec
