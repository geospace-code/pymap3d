""" compute auxiliary latitudes and their inverse"""
from typing import Tuple
import numpy as np
from .ellipsoid import Ellipsoid
from numpy import radians, sin, cos, tan, exp, arctan, log, hypot, degrees, arctan2, sqrt, pi
from .rsphere import rsphere_rectifying
from .rcurve import rcurve_parallel
from .util import wrapToPi, sph2cart, cart2sph

__all__ = ['geodetic2isometric', 'isometric2geodetic', 'geodetic2rectifying', 'rectifying2geodetic', 'geodetic2conformal', 'conformal2geodetic', 'geodetic2parametric', 'parametric2geodetic', 'geodetic2geocentric', 'geocentric2geodetic', 'geodetic2authalic', 'authalic2geodetic']

def geodetic2isometric(lat: float, ell: Ellipsoid = None, deg: bool = True):
    """
     computes isometric latitude of a point on an ellipsoid

     Parameters
     ----------

     lat : float or numpy.ndarray of float
         geodetic latitude
     ell : Ellipsoid, optional
         reference ellipsoid (default WGS84)
     deg : bool, optional
         degrees input/output  (False: radians in/out)

     Returns
     -------

     isolat : float or numpy.ndarray of float
         isometric latiude

     Notes
     -----

     Isometric latitude is an auxiliary latitude proportional to the spacing
     of parallels of latitude on an ellipsoidal mercator projection.

     Based on Deakin, R.E., 2010, 'The Loxodrome on an Ellipsoid', Lecture Notes,
     School of Mathematical and Geospatial Sciences, RMIT University,
     January 2010
    """

    if ell is None:
        ell = Ellipsoid()

    if deg is True:
        lat = np.deg2rad(lat)

    x = ell.eccentricity * np.sin(lat)
    y = (1 - x)/(1 + x)
    z = np.pi/4 + lat/2

#   calculate the isometric latitude
    isolat = np.log(np.tan(z) * (y**(ell.eccentricity /2)))

    if deg is True:
        isolat = np.degrees(isolat)

    return isolat

def isometric2geodetic(isolat, ell: Ellipsoid = None, deg: bool = True):
    '''
    converts from isometric latitude to geodetic latitude
    
    Parameters
     ----------

     isolat : float or numpy.ndarray of float
         isometric latitude
     ell : Ellipsoid, optional
         reference ellipsoid (default WGS84)
     deg : bool, optional
         degrees input/output  (False: radians in/out)

     Returns
     -------

     lat : float or numpy.ndarray of float
         geodetic latiude

     Notes
     -----

    Equations from J. P. Snyder, "Map Projections - A Working Manual",
    US Geological Survey Professional Paper 1395, US Government Printing
    Office, Washington, DC, 1987, pp. 13-18.
'''



    if ell is None:
        ell = Ellipsoid()

    n = ell.thirdflattening
    f1 = 3*n/2 - 27*n**3/32
    f2 = 21*n**2/16 - 55*n**4/32
    f3 = 151*n**3/96
    f4 = 1097*n**4/512

    if deg is True:
        isolat = radians(isolat)

    cnflat = 2 * arctan(exp(isolat)) - (pi/2)
    lat = conformal2geodetic(cnflat, ell, deg=False)

    if deg is True:
        lat = degrees(lat)

    return lat

def geodetic2rectifying(lat, ell: Ellipsoid = None, deg: bool = True):
    '''
    converts from geodetic latitude to rectifying latitude
    
    Parameters
     ----------

     lat : float or numpy.ndarray of float
         geodetic latitude
     ell : Ellipsoid, optional
         reference ellipsoid (default WGS84)
     deg : bool, optional
         degrees input/output  (False: radians in/out)

     Returns
     -------

     reclat : float or numpy.ndarray of float
         rectifying latiude

     Notes
     -----

    Rectifying latitude is an auxiliary latitude used to map an ellipsoid to a 
    sphere such that correct distances along meridians are preserved.

    Equations from J. P. Snyder, "Map Projections - A Working Manual",
    US Geological Survey Professional Paper 1395, US Government Printing
    Office, Washington, DC, 1987, pp. 13-18.
'''

    if ell is None:
        ell = Ellipsoid()

    n = ell.thirdflattening
    f1 = 3*n/2 - 9*n**3/16
    f2 = 15*n**2/16 - 15*n**4/32
    f3 = 35*n**3/48
    f4 = 315*n**4/512

    if deg is True:
        lat = radians(lat)

    reclat = lat - f1*sin(2*lat) + f2*sin(4*lat) - f3*sin(6*lat) + f4*sin(8*lat)

    if deg is True:
        reclat = degrees(reclat)
    
    return reclat


def rectifying2geodetic(reclat, ell: Ellipsoid = None, deg: bool = True):
    '''
    converts from rectifying latitude to geodetic latitude
    
    Parameters
     ----------

     reclat : float or numpy.ndarray of float
         rectifying latitude
     ell : Ellipsoid, optional
         reference ellipsoid (default WGS84)
     deg : bool, optional
         degrees input/output  (False: radians in/out)

     Returns
     -------

     lat : float or numpy.ndarray of float
         geodetic latiude

     Notes
     -----

    Equations from J. P. Snyder, "Map Projections - A Working Manual",
    US Geological Survey Professional Paper 1395, US Government Printing
    Office, Washington, DC, 1987, pp. 13-18.
'''

    if ell is None:
        ell = Ellipsoid()

    n = ell.thirdflattening
    f1 = 3*n/2 - 27*n**3/32
    f2 = 21*n**2/16 - 55*n**4/32
    f3 = 151*n**3/96
    f4 = 1097*n**4/512

    if deg is True:
        reclat = radians(reclat)

    lat = reclat + f1*sin(2*reclat) + f2*sin(4*reclat) + f3*sin(6*reclat) + f4*sin(8*reclat)

    if deg is True:
        lat = degrees(lat)

    return lat


def geodetic2conformal(lat, ell: Ellipsoid = None, deg: bool = True):
    '''
    converts from geodetic latitude to conformal latitude
    
    Parameters
     ----------

     lat : float or numpy.ndarray of float
         geodetic latitude
     ell : Ellipsoid, optional
         reference ellipsoid (default WGS84)
     deg : bool, optional
         degrees input/output  (False: radians in/out)

     Returns
     -------

     cnflat : float or numpy.ndarray of float
         conformal latiude

     Notes
     -----

    Equations from J. P. Snyder, "Map Projections - A Working Manual",
    US Geological Survey Professional Paper 1395, US Government Printing
    Office, Washington, DC, 1987, pp. 13-18.
    '''

    if ell is None:
        ell = Ellipsoid()

    if deg is True:
        lat = radians(lat)

    e = ell.eccentricity
    f1 = 1 - e*sin(lat)
    f2 = 1 + e*sin(lat)
    f3 = 1 - sin(lat)
    f4 = 1 + sin(lat)

    #  compute conformal latitudes with correction for points at +90
    cnflat = np.where(f3 == 0,  pi/2, 2 * arctan(sqrt((f4/f3) * ((f1/f2)**e))) - (pi/2))

    if deg is True:
        cnflat = degrees(cnflat)

    return cnflat

def conformal2geodetic(cnflat, ell: Ellipsoid = None, deg: bool = True):
    '''
    converts from conformal latitude to geodetic latitude

    Parameters
    ----------

    cnflat : float or numpy.ndarray of float
        conformal latitude
    ell : Ellipsoid, optional
        reference ellipsoid (default WGS84)
    deg : bool, optional
        degrees input/output  (False: radians in/out)

    Returns
    -------

    lat : float or numpy.ndarray of float
        geodetic latiude

    Notes
    -----

    Equations from J. P. Snyder, "Map Projections - A Working Manual",
    US Geological Survey Professional Paper 1395, US Government Printing
    Office, Washington, DC, 1987, pp. 13-18.
    '''

    if ell is None:
        ell = Ellipsoid()

    if deg is True:
        cnflat = radians(cnflat)

    e = ell.eccentricity
    f1 = e**2/2 + 5*e**4/24 + e**6/12 + 13*e**8/360
    f2 = 7*e**4/48 + 29*e**6/240 + 811*e**8/11520
    f3 = 7*e**6/120 + 81*e**8/1120
    f4 = 4279*e**8/161280

    lat = cnflat + f1*sin(2*cnflat) + f2*sin(4*cnflat) + f3*sin(6*cnflat) + f4*sin(8*cnflat)

    if deg is True:
        lat = degrees(lat)

    return lat

def geodetic2parametric(lat, ell: Ellipsoid = None, deg: bool = True):
    '''
    converts from geodetic latitude to parametric latitude
    
    Parameters
     ----------

     lat : float or numpy.ndarray of float
         geodetic latitude
     ell : Ellipsoid, optional
         reference ellipsoid (default WGS84)
     deg : bool, optional
         degrees input/output  (False: radians in/out)

     Returns
     -------

     parlat : float or numpy.ndarray of float
         parametric latiude

     Notes
     -----

    Equations from J. P. Snyder, "Map Projections - A Working Manual",
    US Geological Survey Professional Paper 1395, US Government Printing
    Office, Washington, DC, 1987, pp. 13-18.
    '''

    if ell is None:
        ell = Ellipsoid()

    if deg is True:
        lat = radians(lat)

    epsilon = 1E-10

    lat = np.where(lat == pi/2, pi/2 - epsilon, lat)
    lat = np.where(lat == -pi/2, -pi/2 + epsilon, lat)

    parlat = arctan(sqrt(1-(ell.eccentricity)**2) * tan(lat))

    if deg is True:
        parlat = degrees(parlat)

    return parlat

def parametric2geodetic(parlat, ell: Ellipsoid = None, deg: bool = True):
    '''
    converts from parametric latitude to geodetic latitude
    
    Parameters
     ----------

     parlat : float or numpy.ndarray of float
         parametric latitude
     ell : Ellipsoid, optional
         reference ellipsoid (default WGS84)
     deg : bool, optional
         degrees input/output  (False: radians in/out)

     Returns
     -------

     lat : float or numpy.ndarray of float
         geodetic latiude

     Notes
     -----

    Equations from J. P. Snyder, "Map Projections - A Working Manual",
    US Geological Survey Professional Paper 1395, US Government Printing
    Office, Washington, DC, 1987, pp. 13-18.
    '''

    if ell is None:
        ell = Ellipsoid()

    if deg is True:
        parlat = radians(parlat)

    epsilon = 1E-10

    parlat = np.where(parlat == pi/2, pi/2 - epsilon, parlat)
    parlat = np.where(parlat == -pi/2, -pi/2 + epsilon, parlat)

    lat = arctan(tan(parlat)/ sqrt(1-(ell.eccentricity)**2))

    if deg is True:
        lat = degrees(lat)

    return lat

def geodetic2geocentric(lat, ell: Ellipsoid = None, deg: bool = True):
    '''
    converts from geodetic latitude to geocentric latitude
    
    Parameters
     ----------

     lat : float or numpy.ndarray of float
         geodetic latitude
     ell : Ellipsoid, optional
         reference ellipsoid (default WGS84)
     deg : bool, optional
         degrees input/output  (False: radians in/out)

     Returns
     -------

     cenlat : float or numpy.ndarray of float
         geocentric latiude

     Notes
     -----

    Equations from J. P. Snyder, "Map Projections - A Working Manual",
    US Geological Survey Professional Paper 1395, US Government Printing
    Office, Washington, DC, 1987, pp. 13-18.
    '''

    if ell is None:
        ell = Ellipsoid()

    if deg is True:
        lat = radians(lat)

    epsilon = 1E-10

    lat = np.where(lat == pi/2, pi/2 - epsilon, lat)
    lat = np.where(lat == -pi/2, -pi/2 + epsilon, lat)

    cenlat = arctan((1-(ell.eccentricity)**2) * tan(lat))

    if deg is True:
        cenlat = degrees(cenlat)

    return cenlat

def geocentric2geodetic(cenlat, ell: Ellipsoid = None, deg: bool = True):
    '''
    converts from geocentric latitude to geodetic latitude
    
    Parameters
     ----------

     cenlat : float or numpy.ndarray of float
         geocentric latitude
     ell : Ellipsoid, optional
         reference ellipsoid (default WGS84)
     deg : bool, optional
         degrees input/output  (False: radians in/out)

     Returns
     -------

     lat : float or numpy.ndarray of float
         geodetic latiude

     Notes
     -----

    Equations from J. P. Snyder, "Map Projections - A Working Manual",
    US Geological Survey Professional Paper 1395, US Government Printing
    Office, Washington, DC, 1987, pp. 13-18.
    '''

    if ell is None:
        ell = Ellipsoid()

    if deg is True:
        cenlat = radians(cenlat)

    epsilon = 1E-10

    cenlat = np.where(cenlat == pi/2, pi/2 - epsilon, cenlat)
    cenlat = np.where(cenlat == -pi/2, -pi/2 + epsilon, cenlat)

    lat = arctan(tan(cenlat)/ (1-(ell.eccentricity)**2))

    if deg is True:
        lat = degrees(lat)

    return lat

def geodetic2authalic(lat, ell: Ellipsoid = None, deg: bool = True):
    '''
    converts from geodetic latitude to authalic latitude
    
    Parameters
     ----------

     lat : float or numpy.ndarray of float
         geodetic latitude
     ell : Ellipsoid, optional
         reference ellipsoid (default WGS84)
     deg : bool, optional
         degrees input/output  (False: radians in/out)

     Returns
     -------

     autlat : float or numpy.ndarray of float
         authalic latiude

     Notes
     -----

    Equations from J. P. Snyder, "Map Projections - A Working Manual",
    US Geological Survey Professional Paper 1395, US Government Printing
    Office, Washington, DC, 1987, pp. 13-18.
'''

    if ell is None:
        ell = Ellipsoid()

    if deg is True:
        lat = radians(lat)

    e = ell.eccentricity
    f1 = e**2/3 + 31*e**4/180 + 59*e**6/560
    f2 = 17*e**4/360 + 61*e**6/1260
    f3 = 383*e**6/45360

    autlat = lat - f1*sin(2*lat) + f2*sin(4*lat) - f3*sin(6*lat) 

    if deg is True:
        autlat = degrees(autlat)

    return autlat

def authalic2geodetic(autlat, ell: Ellipsoid = None, deg: bool = True):
    '''
    converts from authalic latitude to geodetic latitude
    
    Parameters
     ----------

     autlat : float or numpy.ndarray of float
         authalic latitude
     ell : Ellipsoid, optional
         reference ellipsoid (default WGS84)
     deg : bool, optional
         degrees input/output  (False: radians in/out)

     Returns
     -------

     lat : float or numpy.ndarray of float
         geodetic latiude

     Notes
     -----

    Equations from J. P. Snyder, "Map Projections - A Working Manual",
    US Geological Survey Professional Paper 1395, US Government Printing
    Office, Washington, DC, 1987, pp. 13-18.
'''

    if ell is None:
        ell = Ellipsoid()

    if deg is True:
        autlat = radians(autlat)

    e = ell.eccentricity
    f1 = e**2/3 + 31*e**4/180 + 517*e**6/5040
    f2 = 23*e**4/360 + 251*e**6/3780
    f3 = 761*e**6/45360

    lat = autlat + f1*sin(2*autlat) + f2*sin(4*autlat) + f3*sin(6*autlat)   

    if deg is True:
        lat = degrees(lat)

    return lat