from typing import Tuple
from datetime import datetime
import numpy as np
from .ecef import ecef2enu, geodetic2ecef, ecef2geodetic, enu2uvw
from .enu import geodetic2enu, aer2enu, enu2aer
from .eci import eci2ecef, ecef2eci


def ecef2aer(x: float, y: float, z: float,
             lat0: float, lon0: float, h0: float,
             ell=None, deg: bool = True) -> Tuple[float, float, float]:
    """
    `ecef2aer` gives azimuth, elevation and slant range from an Observer to a Point with ECEF coordinates.

    ## Inputs
    
    * x,y,z  [meters] target ECEF location                             [0,Infinity)
    * lat0, lon0 (degrees/radians)  Observer coordinates on ellipsoid  [-90,90],[-180,180]
    * h0     [meters]                observer altitude
    * ell    reference ellipsoid
    * deg    degrees input/output  (False: radians in/out)

    ## Outputs
    
    * azimuth, elevation (degrees/radians)                             [0,360),[0,90]
    * slant range [meters]                                             [0,Infinity)
    """
    xEast, yNorth, zUp = ecef2enu(x, y, z, lat0, lon0, h0, ell, deg=deg)

    return enu2aer(xEast, yNorth, zUp, deg=deg)


def geodetic2aer(lat: float, lon: float, h: float,
                 lat0: float, lon0: float, h0: float,
                 ell=None, deg: bool = True) -> Tuple[float, float, float]:
    """
    `geodetic2aer` gives azimuth, elevation and slant range from an Observer to a Point with geodetic coordinates.


    ## InputS

    * Target:   lat, lon, h (altitude, meters)
    * Observer: lat0, lon0, h0 (altitude, meters)
    * ell    reference ellipsoid
    * deg    degrees input/output  (False: radians in/out)

    ## Outputs
    
    * azimuth, elevation (degrees/radians)
    * slant range [meters]
    """
    e, n, u = geodetic2enu(lat, lon, h, lat0, lon0, h0, ell, deg=deg)

    return enu2aer(e, n, u, deg=deg)


def aer2geodetic(az: float, el: float, srange: float,
                 lat0: float, lon0: float, h0: float,
                 deg: bool = True) -> Tuple[float, float, float]:
    """
    `aer2geodetic` gives geodetic coordinates of a point with az, el, range from an observer
    
    ## Inputs
    
    * az,el (degrees/radians)
    * srange[meters]        [0, Infinity)

    * Observer: lat0,lon0 [degrees] altitude h0 [meters]

    * deg :   degrees input/output  (False: radians in/out)

    ## Outputs
    * WGS84 lat,lon [degrees]  h0altitude above spheroid  [meters]

    """
    x, y, z = aer2ecef(az, el, srange, lat0, lon0, h0, deg=deg)

    return ecef2geodetic(x, y, z, deg=deg)


def eci2aer(eci: Tuple[float, float, float],
            lat0: float, lon0: float, h0: float,
            t: datetime,
            useastropy: bool = True) -> Tuple[float, float, float]:
    """
    `eci2aer` takes ECI coordinates of point and gives az, el, slant range from Observer

    ## Inputs

    * eci [meters] Nx3 target ECI location (x,y,z)                    [0,Infinity)
    * lat0, lon0 (degrees/radians)  Observer coordinates on ellipsoid [-90,90],[-180,180]
    * h0     [meters]                observer altitude                [0,Infinity)
    * t  time (datetime.datetime)   time of obsevation (UTC)


    ## Outputs
    
    * azimuth, elevation (degrees/radians)                             [0,360),[0,90]
    * slant range [meters]                                             [0,Infinity)
    """
    ecef = np.atleast_2d(eci2ecef(eci, t, useastropy))

    return ecef2aer(ecef[:, 0], ecef[:, 1], ecef[:, 2], lat0, lon0, h0)


def aer2eci(az: float, el: float, srange: float,
            lat0: float, lon0: float, h0: float, t: datetime,
            ell=None, deg: bool = True,
            useastropy: bool = True) -> np.ndarray:
    """
    `aer2eci` gives ECI of a point from an observer at az, el, slant range
    
    ## Inputs
    
    * azimuth, elevation (degrees/radians)                             [0,360),[0,90]
    * slant range [meters]                                             [0,Infinity)
    * Observer: lat0, lon0, h0 (altitude, meters)
    * ell    reference ellipsoid
    * deg    degrees input/output  (False: radians in/out)
    * t  datetime.datetime of obseration

    ## Outputs
    
    * eci  x,y,z (meters)
    """
    x, y, z = aer2ecef(az, el, srange, lat0, lon0, h0, ell, deg)

    return ecef2eci(np.column_stack((x, y, z)), t, useastropy)


def aer2ecef(az: float, el: float, srange: float,
             lat0: float, lon0: float, alt0: float,
             ell=None, deg: bool = True) -> Tuple[float, float, float]:
    """
    `aer2ecef` converts target azimuth, elevation, range (meters) from observer at lat0,lon0,alt0 to ECEF coordinates.

    ## Inputs
     
    * azimuth, elevation (degrees/radians)                             [0,360),[0,90]
    * slant range [meters]                                             [0,Infinity)
    * Observer: lat0, lon0, h0 (altitude, meters)
    * ell    reference ellipsoid
    * deg    degrees input/output  (False: radians in/out)

    ## OutputS
    
    ECEF x,y,z  [meters]

    if you specify NaN for srange, return value z will be NaN
    """
    # Origin of the local system in geocentric coordinates.
    x0, y0, z0 = geodetic2ecef(lat0, lon0, alt0, ell, deg=deg)
    # Convert Local Spherical AER to ENU
    e1, n1, u1 = aer2enu(az, el, srange, deg=deg)
    # Rotating ENU to ECEF
    dx, dy, dz = enu2uvw(e1, n1, u1, lat0, lon0, deg=deg)
    # Origin + offset from origin equals position in ECEF
    return x0 + dx, y0 + dy, z0 + dz
