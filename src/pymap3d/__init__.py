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

Most functions accept scalar float or NumPy NDArray.

Other languages
---------------

Companion packages exist for:

* Matlab / GNU Octave: [Matmap3D](https://github.com/geospace-code/matmap3d)
* Fortran: [Maptran3D](https://github.com/geospace-code/maptran3d)
"""

__version__ = "3.2.0"

from .aer import aer2ecef, aer2geodetic, ecef2aer, geodetic2aer
from .body import (
    body2ecef,
    body2ecefv,
    body2enu,
    body2enu_matrix,
    body2geodetic,
    body2ned,
    body2ned_matrix,
    body2nwu,
    body2nwu_matrix,
    body2sez,
    body2sez_matrix,
    body_matrix,
    ecef2body,
    ecef2bodyv,
    enu2body,
    enu2body_matrix,
    geodetic2body,
    ned2body,
    ned2body_matrix,
    nwu2body,
    nwu2body_matrix,
    sez2body,
    sez2body_matrix,
)
from .dca import (
    enu2dca,
    dca2enu,
    dca2ned,
    ned2dca,
    ecef2dca,
    dca2ecef,
    geodetic2dca,
    dca2geodetic,
    aer2dca,
    dca2aer,
)

from .ecef import (
    ecef2enu,
    ecef2enuv,
    ecef2geodetic,
    eci2geodetic,
    enu2ecef,
    enu2uvw,
    geodetic2ecef,
    geodetic2eci,
    uvw2enu,
)
from .ellipsoid import Ellipsoid
from .enu import (
    aer2enu,
    enu2aer,
    enu2geodetic,
    geodetic2enu,
    enu2ecefv,
    ecef2enu_matrix,
    enu2ecef_matrix,
)
from .ned import (
    aer2ned,
    ecef2ned,
    ecef2nedv,
    ecef2ned_matrix,
    geodetic2ned,
    ned2ecef_matrix,
    ned2aer,
    ned2ecef,
    ned2geodetic,
)
from .nwu import (
    aer2nwu,
    ecef2nwu,
    ecef2nwuv,
    ecef2nwu_matrix,
    geodetic2nwu,
    nwu2aer,
    nwu2ecef,
    nwu2ecefv,
    nwu2ecef_matrix,
    nwu2geodetic,
)
from .sez import (
    aer2sez,
    ecef2sez,
    ecef2sezv,
    ecef2sez_matrix,
    geodetic2sez,
    sez2aer,
    sez2ecef,
    sez2ecefv,
    sez2ecef_matrix,
    sez2geodetic,
)
from .sidereal import datetime2sidereal, greenwichsrt
from .spherical import geodetic2spherical, spherical2geodetic
from .timeconv import str2dt

from .latitude import (
    geodetic2isometric,
    isometric2geodetic,
    geodetic2rectifying,
    rectifying2geodetic,
    geodetic2conformal,
    conformal2geodetic,
    geodetic2parametric,
    parametric2geodetic,
    geodetic2geocentric,
    geocentric2geodetic,
    geodetic2authalic,
    authalic2geodetic,
    geod2geoc,
    geoc2geod,
)

from .nvector import (
    geodetic2nvector,
    nvector2geodetic,
    ecef2nvector,
    nvector2ecef,
    nvector_distance,
    nvector_interpolate,
    nvector_mean,
    nvector_cross_track_distance,
    nvector_intersection,
)

from .rcurve import parallel, meridian, transverse, geocentric_radius

__all__ = [
    "aer2ecef",
    "aer2geodetic",
    "ecef2aer",
    "geodetic2aer",
    "body_matrix",
    "ned2body_matrix",
    "body2ned_matrix",
    "enu2body_matrix",
    "body2enu_matrix",
    "ned2body",
    "body2ned",
    "enu2body",
    "body2enu",
    "sez2body",
    "body2sez",
    "nwu2body",
    "body2nwu",
    "sez2body_matrix",
    "body2sez_matrix",
    "nwu2body_matrix",
    "body2nwu_matrix",
    "ecef2body",
    "body2ecef",
    "ecef2bodyv",
    "body2ecefv",
    "geodetic2body",
    "body2geodetic",
    "ecef2enu",
    "ecef2enuv",
    "ecef2geodetic",
    "eci2geodetic",
    "enu2ecef",
    "enu2uvw",
    "geodetic2ecef",
    "geodetic2eci",
    "uvw2enu",
    "Ellipsoid",
    "aer2enu",
    "enu2aer",
    "enu2geodetic",
    "enu2ecefv",
    "geodetic2enu",
    "ecef2enu_matrix",
    "enu2ecef_matrix",
    "aer2ned",
    "ecef2ned",
    "ecef2nedv",
    "ecef2ned_matrix",
    "geodetic2ned",
    "ned2aer",
    "ned2ecef",
    "ned2ecef_matrix",
    "ned2geodetic",
    "aer2sez",
    "ecef2sez",
    "ecef2sezv",
    "ecef2sez_matrix",
    "geodetic2sez",
    "sez2aer",
    "sez2ecef",
    "sez2ecefv",
    "sez2ecef_matrix",
    "sez2geodetic",
    "aer2nwu",
    "ecef2nwu",
    "ecef2nwuv",
    "ecef2nwu_matrix",
    "geodetic2nwu",
    "nwu2aer",
    "nwu2ecef",
    "nwu2ecefv",
    "nwu2ecef_matrix",
    "nwu2geodetic",
    "datetime2sidereal",
    "greenwichsrt",
    "geodetic2spherical",
    "spherical2geodetic",
    "str2dt",
    "azel2radec",
    "radec2azel",
    "parallel",
    "meridian",
    "transverse",
    "geocentric_radius",
    "geodetic2isometric",
    "isometric2geodetic",
    "geodetic2rectifying",
    "rectifying2geodetic",
    "geodetic2conformal",
    "conformal2geodetic",
    "geodetic2parametric",
    "parametric2geodetic",
    "geodetic2geocentric",
    "geocentric2geodetic",
    "geodetic2authalic",
    "authalic2geodetic",
    "geod2geoc",
    "geoc2geod",
    "nvector2geodetic",
    "geodetic2nvector",
    "ecef2nvector",
    "nvector2ecef",
    "enu2dca",
    "dca2enu",
    "dca2ned",
    "ned2dca",
    "ecef2dca",
    "dca2ecef",
    "geodetic2dca",
    "dca2geodetic",
    "aer2dca",
    "dca2aer",
    "nvector_distance",
    "nvector_interpolate",
    "nvector_mean",
    "nvector_cross_track_distance",
    "nvector_intersection",
]


from .aer import aer2eci, eci2aer
from .azelradec import azel2radec, radec2azel
from .eci import ecef2eci, ecef2eci_state, eci2ecef, eci2ecef_state

__all__ += [
    "aer2eci",
    "eci2aer",
    "ecef2eci",
    "eci2ecef",
    "ecef2eci_state",
    "eci2ecef_state",
]
