""" isometric latitude, meridian distance """
import numpy as np
from .ecef import Ellipsoid


def isometric(lat: float, ell: Ellipsoid = None, deg: bool = True):
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

    f = ell.f  # flattening of ellipsoid

    if deg is True:
        lat = np.deg2rad(lat)

    e2 = f * (2 - f)  # eccentricity-squared
    e = np.sqrt(e2)  # eccentricity of ellipsoid

    x = e * np.sin(lat)
    y = (1 - x) / (1 + x)
    z = np.pi / 4 + lat / 2

#   calculate the isometric latitude
    isolat = np.log(np.tan(z) * (y**(e / 2)))

    if deg is True:
        isolat = np.degrees(isolat)

    return isolat


def meridian_dist(lat: float, ell: Ellipsoid = None, deg: bool = True):
    """
    computes the meridian distance on an ellipsoid *from the equator* to a latitude.

    Parameters
    ----------

    lat : float or numpy.ndarray of float
        geodetic latitude
    ell : Ellipsoid, optional
         reference ellipsoid (default WGS84)
    deg : bool, optional
         degrees input/output  (False: radians in/out)

    Results
    -------

    mdist : float or numpy.ndarray of float
         meridian distance (degrees/radians)


    Notes
    -----

    Formula given Baeschlin, C.F., 1948,
    "Lehrbuch Der Geodasie", Orell Fussli Verlag, Zurich, pp.47-50.

    Based on Deakin, R.E., 2010, 'The Loxodrome on an Ellipsoid', Lecture Notes,
    School of Mathematical and Geospatial Sciences, RMIT University, January 2010
    """

    if deg is True:
        lat = np.radians(lat)

    #   set ellipsoid parameters
    if ell is None:
        ell = Ellipsoid()

    a = ell.a
    f = ell.f  # flattening of ellipsoid

    e2 = f * (2 - f)  # eccentricity-squared

    # powers of eccentricity
    e4 = e2 * e2
    e6 = e4 * e2
    e8 = e6 * e2
    e10 = e8 * e2

    # coefficients of series expansion for meridian distance
    A = 1 + (3 / 4) * e2 + (45 / 64) * e4 + (175 / 256) * e6 + (11025 / 16384) * e8 + (43659 / 65536) * e10
    B = (3 / 4) * e2 + (15 / 16) * e4 + (525 / 512) * e6 + (2205 / 2048) * e8 + (72765 / 65536) * e10
    C = (15 / 64) * e4 + (105 / 256) * e6 + (2205 / 4096) * e8 + (10395 / 16384) * e10
    D = (35 / 512) * e6 + (315 / 2048) * e8 + (31185 / 131072) * e10
    E = (315 / 16384) * e8 + (3465 / 65536) * e10
    F = (693 / 131072) * e10

    term1 = A * lat
    term2 = (B / 2) * np.sin(2 * lat)
    term3 = (C / 4) * np.sin(4 * lat)
    term4 = (D / 6) * np.sin(6 * lat)
    term5 = (E / 8) * np.sin(8 * lat)
    term6 = (F / 10) * np.sin(10 * lat)

    mdist = a * (1 - e2) * (term1 - term2 + term3 - term4 + term5 - term6)

    return mdist


def loxodrome_inverse(lat1: float, lon1: float, lat2: float, lon2: float,
                      ell: Ellipsoid = None, deg: bool = True):
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

    lox_s : float or numpy.ndarray of float
        distance along loxodrome
    az12 : float or numpy.ndarray of float
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
    isolat1 = isometric(lat1, deg=False, ell=ell)
    isolat2 = isometric(lat2, deg=False, ell=ell)

    # compute changes in isometric latitude and longitude between points
    disolat = isolat2 - isolat1
    dlon = lon2 - lon1

    # compute azimuth
    az12 = np.arctan2(dlon, disolat)

    # compute distance along loxodromic curve
    m1 = meridian_dist(lat1, deg=False, ell=ell)
    m2 = meridian_dist(lat2, deg=False, ell=ell)
    dm = m2 - m1
    lox_s = dm / np.cos(az12)

    if deg is True:
        az12 = np.degrees(az12) % 360.

    return lox_s, az12
