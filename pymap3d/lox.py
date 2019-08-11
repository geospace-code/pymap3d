import numpy as np
from .ellipsoid import Ellipsoid
from typing import Tuple
from .latitude import geodetic2isometric, geodetic2rectifying, rectifying2geodetic, geodetic2authalic, authalic2geodetic
from .rsphere import rsphere_rectifying
from .rcurve import rcurve_parallel
from .util import wrapToPi, sph2cart, cart2sph
from numpy import arctan, sqrt, tan, sign, sin, cos, arctan2, arcsin, nan, pi, radians, degrees

__all__ = ['loxodrome_inverse', 'loxodrome_direct', 'medianarc', 'departure', 'meanm']

def loxodrome_inverse(lat1: float, lon1: float, lat2: float, lon2: float,
                      ell: Ellipsoid = None, deg: bool = True)  -> Tuple[float, float]:
    """
    computes the arc length and azimuth of the loxodrome
    between two points on the surface of the reference ellipsoid

    Parameters
    ----------

    lat1 : float or numpy.ndarray of float
        geodetic latitude of first point
    lon1 : float or numpy.ndarray of float
        geodetic longitude of first point
    lat2 : float or numpy.ndarray of float
        geodetic latitude of second point
    lon2 : float or numpy.ndarray of float
        geodetic longitude of second point
    ell : Ellipsoid, optional
         reference ellipsoid (default WGS84)
    deg : bool, optional
         degrees input/output  (False: radians in/out)

    Results
    -------

    dist : float or numpy.ndarray of float
        distance along loxodrome
    a12 : float or numpy.ndarray of float
        azimuth of loxodrome (degrees/radians)

    Based on Deakin, R.E., 2010, 'The Loxodrome on an Ellipsoid', Lecture Notes,
    School of Mathematical and Geospatial Sciences, RMIT University, January 2010

    [1] Bowring, B.R., 1985, 'The geometry of the loxodrome on the
    ellipsoid', The Canadian Surveyor, Vol. 39, No. 3, Autumn 1985,
    pp.223-230.
    [2] Snyder, J.P., 1987, Map Projections-A Working Manual. U.S.
    Geological Survey Professional Paper 1395. Washington, DC: U.S.
    Government Printing Office, pp.15-16 and pp. 44-45.
    [3] Thomas, P.D., 1952, Conformal Projections in Geodesy and
    Cartography, Special Publication No. 251, Coast and Geodetic
    Survey, U.S. Department of Commerce, Washington, DC: U.S.
    Government Printing Office, p. 66.
    """

    #   set ellipsoid parameters
    if ell is None:
        ell = Ellipsoid()

    if deg is True:
        lat1, lon1, lat2, lon2 = np.radians([lat1, lon1, lat2, lon2])

    # compute isometric latitude of P1 and P2
    isolat1 = geodetic2isometric(lat1, deg=False, ell=ell)
    isolat2 = geodetic2isometric(lat2, deg=False, ell=ell)

    # compute changes in isometric latitude and longitude between points
    disolat = isolat2 - isolat1

    dlon = abs(wrapToPi(lon2-lon1))

    # compute azimuth
    a12 = np.arctan2(dlon, disolat)
    cosaz = abs(cos(a12))

    # compute distance along loxodromic curve
    dist = meridianarc(lat2, lat1, deg=False, ell=ell) / cosaz

    # consider degenerate case (directly east/west)
    epsilon = 1e-10 
    dist =  np.where(cosaz < epsilon,departure(lon1, lon2, (lat1+lat2)/2, ell, deg=False), dist)

    if deg is True:
        a12 = np.degrees(a12) % 360.

    return dist, a12


def loxodrome_direct(lat1: float, lon1: float, rng: float, a12: float,
            ell: Ellipsoid = None, deg: bool = True) -> Tuple[float, float, float]:
    """

    Parameters
    ----------

    lat1 : float
        inital geodetic latitude (degrees)
    lon1 : float
        initial geodetic longitude (degrees)
    rng : float
        ground distance (meters)
    a12 : float
        azimuth (degrees) clockwide from north.
    ell : Ellipsoid, optional
          reference ellipsoid

    Results
    -------

    lat2 : float
        final geodetic latitude (degrees)
    lon2 : float
        final geodetic longitude (degrees)
    a21 : float
        reverse azimuth (degrees), at final point facing back toward the intial point

"""
    if ell is None:
        ell = Ellipsoid()

    lat1 = np.atleast_1d(lat1)
    lon1 = np.atleast_1d(lon1)
    rng = np.atleast_1d(rng)
    a12 = np.atleast_1d(a12)

    if rng.ndim != 1 or a12.ndim != 1:
        raise ValueError('Range and azimuth must be scalar or vector')

    if lat1.size > 1 and rng.size > 1:
        raise ValueError('LOXODROME_DIRECT: Variable ranges are only allowed for a single point.')

    if rng.size != a12.size and rng.size == 1:
        rng = np.broadcast_to(rng, a12.size)

    if deg is True:
        lat1, lon1, a12 = radians([lat1, lon1, a12])
    
    if abs(lat1) > (pi/2):
        raise ValueError('LOXODROME_DIRECT: Input lat. must be between -90 and 90 deg., inclusive.')

    # if range is negative, use inverse bearing
    a12 = np.where(sign(rng) == -1, (a12 + pi) % (2*pi), a12)
    rng = abs(rng)

    #   compute rectifying sphere latitude and radius
    reclat = geodetic2rectifying(lat1, ell, deg=False)

    lat2  = np.zeros_like(lat1)    
    lon2  = np.zeros_like(lon1)
    epsilon = 1e-10     # Set tolerance (should set this for specific machine)

    # correct azimuths at either pole.
    a12 = np.where(lat1 >= (pi/2) - epsilon, pi, a12)
    a12 = np.where(lat1 <= (-pi/2) + epsilon, 0,  a12)

    # compute the new points
    cosaz  = cos(a12)
    lat2 = reclat + (rng/rsphere_rectifying(ell))*cosaz # compute rectifying latitude
    lat2 = rectifying2geodetic(lat2,ell,deg=False) #  transform to geodetic latitude

    # nudge latitudes near either pole
    lat2 = np.where(lat2 >=  (pi/2) - epsilon, (pi/2) - epsilon, lat2)
    lat2 = np.where(lat2 <=  (-pi/2) + epsilon, (-pi/2) + epsilon , lat2)
    lat2 = np.where(abs(cosaz)<=epsilon, lat1, lat2)

    newiso = geodetic2isometric(lat2,ell,deg=False)
    iso    = geodetic2isometric(lat1,ell,deg=False)
    
    dlon = np.where(abs(cosaz)<=epsilon, sign(sin(a12)) * rng/rcurve_parallel(lat1, ell, deg=False), tan(a12) * (newiso - iso))
    lon2 = lon1 + dlon

    lon2 = wrapToPi(lon2)
    a21 = (a12 + pi)%(2*pi)

    if deg is True:
        lat2, lon2, a21 = degrees([lat2, lon2, a21])

    return lat2, lon2, a21

def meridianarc(lat1: float, lat2: float, ell: Ellipsoid = None, deg: bool = True) -> float:
    '''
    Computes the meridian distance on an ellipsoid between two latitudes.

    Parameters
    ----------

    lat1, lat2 : float or numpy.ndarray of float
        geodetic latitudes
    ell : Ellipsoid, optional
         reference ellipsoid (default WGS84)
    deg : bool, optional
         degrees input/output  (False: radians in/out)

    Results
    -------

    dist : float or numpy.ndarray of float
         distance (units same as ellipsoid)
    '''

    if deg is True:
        lat1, lat2 = np.radians([lat1, lat2])

    #   set ellipsoid parameters
    if ell is None:
        ell = Ellipsoid()
    
    rlat1 = geodetic2rectifying(lat1, ell, deg=False)
    rlat2 = geodetic2rectifying(lat2, ell, deg=False)

    return rsphere_rectifying(ell) * abs(rlat2 - rlat1)

def departure(lon1: float, lon2: float, lat: float, ell: Ellipsoid = None, deg: bool = True) -> float:
    '''
    Computes is the distance along a specific parallel between two meridians.
    '''

    if ell is None:
        ell = Ellipsoid()

    if deg is True:
        lon1, lon2, lat = np.radians([lon1, lon2, lat])
    
    return rcurve_parallel(lat, ell, deg=False) * abs(wrapToPi(lon2-lon1))


def meanm(lat: float, lon: float, ell: Ellipsoid = None, deg: bool = True) -> Tuple[float, float]:
    '''Computes geographic mean for geographic points on an ellipsoid'''
    
    if ell is None:
        ell = Ellipsoid()

    if deg is True:
        lon, lat = radians([lon, lat])

    lat = geodetic2authalic(lat)
    x,y,z = sph2cart(lon,lat,np.ones_like(lat))
    latbar, lonbar, _ =cart2sph(np.sum(x),np.sum(y),np.sum(z))
    latbar = authalic2geodetic(latbar, ell, deg=False)
    lonbar = wrapToPi(lonbar)

    return latbar, lonbar
