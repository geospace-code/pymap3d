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
import sys

if sys.version_info >= (3, 5):
    try:
        import numpy
        import dateutil  # noqa: F401
    except ImportError:
        numpy = None
else:
    numpy = None

if numpy is not None:
    from .timeconv import str2dt  # noqa: F401
    from .azelradec import radec2azel, azel2radec  # noqa: F401
    from .eci import eci2ecef, ecef2eci  # noqa: F401
    from .sidereal import datetime2sidereal  # noqa: F401
    from .ecef import geodetic2ecef, ecef2geodetic, eci2geodetic, ecef2enuv, enu2ecef, ecef2enu, enu2uvw, uvw2enu  # noqa: F401
    from .ellipsoid import Ellipsoid  # noqa: F401
    from .ned import ned2ecef, ned2geodetic, geodetic2ned, ecef2nedv, ned2aer, aer2ned, ecef2ned  # noqa: F401
    from .enu import enu2geodetic, geodetic2enu, aer2enu, enu2aer  # noqa: F401
    from .aer import ecef2aer, aer2ecef, geodetic2aer, aer2geodetic, eci2aer, aer2eci  # noqa: F401
    from .los import lookAtSpheroid  # noqa: F401
    from .lox import loxodrome_direct, loxodrome_inverse, meridianarc, departure, meanm
    from .latitude import geodetic2isometric, isometric2geodetic, geodetic2rectifying, rectifying2geodetic, geodetic2conformal, conformal2geodetic, geodetic2parametric, parametric2geodetic, geodetic2geocentric, geocentric2geodetic, geodetic2authalic, authalic2geodetic
    from .rcurve import rcurve_meridian, rcurve_parallel, rcurve_transverse
    from .rsphere import rsphere_authalic, rsphere_biaxial, rsphere_curve, rsphere_eqavol, rsphere_euler, rsphere_rectifying, rsphere_triaxial
    from .util import wrapToPi, wrapTo2Pi, sph2cart, cart2sph, pol2cart, cart2pol
else:  # pure Python only
    from .math import *  # type: ignore  # noqa: F401, F403
