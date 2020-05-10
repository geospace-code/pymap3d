""" transforms involving AER: azimuth, elevation, slant range"""
import typing
from datetime import datetime

from .ecef import ecef2enu, geodetic2ecef, ecef2geodetic, enu2uvw
from .enu import geodetic2enu, aer2enu, enu2aer
from .ellipsoid import Ellipsoid

try:
    from .eci import eci2ecef, ecef2eci
except ImportError:
    eci2ecef = ecef2eci = None

__all__ = ["aer2ecef", "ecef2aer", "geodetic2aer", "aer2geodetic", "eci2aer", "aer2eci"]

if typing.TYPE_CHECKING:
    from numpy import ndarray


def ecef2aer(
    x: "ndarray",
    y: "ndarray",
    z: "ndarray",
    lat0: "ndarray",
    lon0: "ndarray",
    h0: "ndarray",
    ell: Ellipsoid = None,
    deg: bool = True,
) -> typing.Tuple["ndarray", "ndarray", "ndarray"]:
    """
    compute azimuth, elevation and slant range from an Observer to a Point with ECEF coordinates.

    ECEF input location is with units of meters

    Parameters
    ----------

    x : "ndarray"
        ECEF x coordinate (meters)
    y : "ndarray"
        ECEF y coordinate (meters)
    z : "ndarray"
        ECEF z coordinate (meters)
    lat0 : "ndarray"
        Observer geodetic latitude
    lon0 : "ndarray"
        Observer geodetic longitude
    h0 : "ndarray"
         observer altitude above geodetic ellipsoid (meters)
    ell : Ellipsoid, optional
        reference ellipsoid
    deg : bool, optional
        degrees input/output  (False: radians in/out)

    Returns
    -------
    az : "ndarray"
         azimuth to target
    el : "ndarray"
         elevation to target
    srange : "ndarray"
         slant range [meters]
    """
    xEast, yNorth, zUp = ecef2enu(x, y, z, lat0, lon0, h0, ell, deg=deg)

    return enu2aer(xEast, yNorth, zUp, deg=deg)


def geodetic2aer(
    lat: "ndarray",
    lon: "ndarray",
    h: "ndarray",
    lat0: "ndarray",
    lon0: "ndarray",
    h0: "ndarray",
    ell: Ellipsoid = None,
    deg: bool = True,
) -> typing.Tuple["ndarray", "ndarray", "ndarray"]:
    """
    gives azimuth, elevation and slant range from an Observer to a Point with geodetic coordinates.


    Parameters
    ----------

    lat : "ndarray"
        target geodetic latitude
    lon : "ndarray"
        target geodetic longitude
    h : "ndarray"
        target altitude above geodetic ellipsoid (meters)
    lat0 : "ndarray"
        Observer geodetic latitude
    lon0 : "ndarray"
        Observer geodetic longitude
    h0 : "ndarray"
         observer altitude above geodetic ellipsoid (meters)
    ell : Ellipsoid, optional
          reference ellipsoid
    deg : bool, optional
          degrees input/output  (False: radians in/out)

    Returns
    -------
    az : "ndarray"
         azimuth
    el : "ndarray"
         elevation
    srange : "ndarray"
         slant range [meters]
    """
    e, n, u = geodetic2enu(lat, lon, h, lat0, lon0, h0, ell, deg=deg)

    return enu2aer(e, n, u, deg=deg)


def aer2geodetic(
    az: "ndarray",
    el: "ndarray",
    srange: "ndarray",
    lat0: "ndarray",
    lon0: "ndarray",
    h0: "ndarray",
    ell: Ellipsoid = None,
    deg: bool = True,
) -> typing.Tuple["ndarray", "ndarray", "ndarray"]:
    """
    gives geodetic coordinates of a point with az, el, range
    from an observer at lat0, lon0, h0

    Parameters
    ----------
    az : "ndarray"
         azimuth to target
    el : "ndarray"
         elevation to target
    srange : "ndarray"
         slant range [meters]
    lat0 : "ndarray"
           Observer geodetic latitude
    lon0 : "ndarray"
           Observer geodetic longitude
    h0 : "ndarray"
         observer altitude above geodetic ellipsoid (meters)
    ell : Ellipsoid, optional
          reference ellipsoid
    deg : bool, optional
          degrees input/output  (False: radians in/out)

    Returns
    -------

    In reference ellipsoid system:

    lat : "ndarray"
          geodetic latitude
    lon : "ndarray"
          geodetic longitude
    alt : "ndarray"
          altitude above ellipsoid  (meters)
    """
    x, y, z = aer2ecef(az, el, srange, lat0, lon0, h0, ell=ell, deg=deg)

    return ecef2geodetic(x, y, z, ell=ell, deg=deg)


def eci2aer(
    x: "ndarray",
    y: "ndarray",
    z: "ndarray",
    lat0: "ndarray",
    lon0: "ndarray",
    h0: "ndarray",
    t: datetime,
    *,
    deg: bool = True,
    use_astropy: bool = True
) -> typing.Tuple["ndarray", "ndarray", "ndarray"]:
    """
    takes Earth Centered Inertial x,y,z ECI coordinates of point and gives az, el, slant range from Observer

    Parameters
    ----------

    x : "ndarray"
        ECI x-location [meters]
    y : "ndarray"
        ECI y-location [meters]
    z : "ndarray"
        ECI z-location [meters]
    lat0 : "ndarray"
           Observer geodetic latitude
    lon0 : "ndarray"
           Observer geodetic longitude
    h0 : "ndarray"
         observer altitude above geodetic ellipsoid (meters)
    t : datetime.datetime
        Observation time
    deg : bool, optional
        true: degrees, false: radians
    use_astropy: bool, optional
        use Astropy (recommended)

    Returns
    -------
    az : "ndarray"
         azimuth to target
    el : "ndarray"
         elevation to target
    srange : "ndarray"
         slant range [meters]
    """
    if eci2ecef is None:
        raise ImportError("pip install numpy")

    xecef, yecef, zecef = eci2ecef(x, y, z, t, use_astropy=use_astropy)

    return ecef2aer(xecef, yecef, zecef, lat0, lon0, h0, deg=deg)


def aer2eci(
    az: "ndarray",
    el: "ndarray",
    srange: "ndarray",
    lat0: "ndarray",
    lon0: "ndarray",
    h0: "ndarray",
    t: datetime,
    ell=None,
    *,
    deg: bool = True,
    use_astropy: bool = True
) -> typing.Tuple["ndarray", "ndarray", "ndarray"]:
    """
    gives ECI of a point from an observer at az, el, slant range

    Parameters
    ----------
    az : "ndarray"
         azimuth to target
    el : "ndarray"
         elevation to target
    srange : "ndarray"
         slant range [meters]
    lat0 : "ndarray"
           Observer geodetic latitude
    lon0 : "ndarray"
           Observer geodetic longitude
    h0 : "ndarray"
         observer altitude above geodetic ellipsoid (meters)
    t : datetime.datetime
        Observation time
    ell : Ellipsoid, optional
          reference ellipsoid
    deg : bool, optional
          degrees input/output  (False: radians in/out)
    use_astropy : bool, optional
        use AstroPy (recommended)

    Returns
    -------

    Earth Centered Inertial x,y,z

    x : "ndarray"
        ECEF x coordinate (meters)
    y : "ndarray"
        ECEF y coordinate (meters)
    z : "ndarray"
        ECEF z coordinate (meters)
    """
    if ecef2eci is None:
        raise ImportError("pip install numpy")

    x, y, z = aer2ecef(az, el, srange, lat0, lon0, h0, ell, deg=deg)

    return ecef2eci(x, y, z, t, use_astropy=use_astropy)


def aer2ecef(
    az: "ndarray",
    el: "ndarray",
    srange: "ndarray",
    lat0: "ndarray",
    lon0: "ndarray",
    alt0: "ndarray",
    ell: Ellipsoid = None,
    deg: bool = True,
) -> typing.Tuple["ndarray", "ndarray", "ndarray"]:
    """
    converts target azimuth, elevation, range from observer at lat0,lon0,alt0 to ECEF coordinates.

    Parameters
    ----------
    az : "ndarray"
         azimuth to target
    el : "ndarray"
         elevation to target
    srange : "ndarray"
         slant range [meters]
    lat0 : "ndarray"
           Observer geodetic latitude
    lon0 : "ndarray"
           Observer geodetic longitude
    h0 : "ndarray"
         observer altitude above geodetic ellipsoid (meters)
    ell : Ellipsoid, optional
          reference ellipsoid
    deg : bool, optional
          degrees input/output  (False: radians in/out)

    Returns
    -------

    ECEF (Earth centered, Earth fixed)  x,y,z

    x : "ndarray"
        ECEF x coordinate (meters)
    y : "ndarray"
        ECEF y coordinate (meters)
    z : "ndarray"
        ECEF z coordinate (meters)


    Notes
    ------
    if srange==NaN, z=NaN
    """
    # Origin of the local system in geocentric coordinates.
    x0, y0, z0 = geodetic2ecef(lat0, lon0, alt0, ell, deg=deg)
    # Convert Local Spherical AER to ENU
    e1, n1, u1 = aer2enu(az, el, srange, deg=deg)
    # Rotating ENU to ECEF
    dx, dy, dz = enu2uvw(e1, n1, u1, lat0, lon0, deg=deg)
    # Origin + offset from origin equals position in ECEF
    return x0 + dx, y0 + dy, z0 + dz
