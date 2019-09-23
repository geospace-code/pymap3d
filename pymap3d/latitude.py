from .ellipsoid import Ellipsoid
from math import atan, radians, degrees, tan

try:
    import numpy
except ImportError:
    numpy = None


def geodetic2geocentric(geodetic_lat: float, ell: Ellipsoid = None, deg: bool = True) -> float:
    if numpy is not None:
        fun = numpy.vectorize(geodetic2geocentric_point)
        return fun(geodetic_lat, ell, deg)
    else:
        return geodetic2geocentric_point(geodetic_lat, ell, deg)


def geodetic2geocentric_point(geodetic_lat: float, ell: Ellipsoid = None, deg: bool = True) -> float:
    """
    convert geodetic latitude to geocentric latitude.

    this is like Matlab geocentricLatitude()
    https://www.mathworks.com/help/map/ref/geocentriclatitude.html

    Parameters
    ----------
    geodetic_lat : float
        geodetic latitude
    ell : Ellipsoid, optional
         reference ellipsoid (default WGS84)
    deg : bool, optional
         degrees input/output  (False: radians in/out)

    Returns
    -------
    geocentric_lat : float
         geocentric latiude

    Notes
    -----
    Equations from J. P. Snyder, "Map Projections - A Working Manual",
    US Geological Survey Professional Paper 1395, US Government Printing
    Office, Washington, DC, 1987, pp. 13-18.
    """

    if ell is None:
        ell = Ellipsoid()

    if abs(geodetic_lat) > 90:
        raise ValueError("-90 <= latitude <= 90")

    if deg is True:
        geodetic_lat = radians(geodetic_lat)

    geocentric_lat = atan((1 - (ell.eccentricity) ** 2) * tan(geodetic_lat))

    if deg is True:
        geocentric_lat = degrees(geocentric_lat)

    return geocentric_lat


def geocentric2geodetic(geocentric_lat: float, ell: Ellipsoid = None, deg: bool = True) -> float:
    if numpy is not None:
        fun = numpy.vectorize(geocentric2geodetic_point)
        return fun(geocentric_lat, ell, deg)
    else:
        return geocentric2geodetic_point(geocentric_lat, ell, deg)


def geocentric2geodetic_point(geocentric_lat: float, ell: Ellipsoid = None, deg: bool = True) -> float:
    """
    converts from geocentric latitude to geodetic latitude

    like Matlab geodeticLatitudeFromGeocentric
    https://www.mathworks.com/help/map/ref/geodeticlatitudefromgeocentric.html

    Parameters
    ----------
    geocentric_lat : float or numpy.ndarray of float
         geocentric latitude
    ell : Ellipsoid, optional
         reference ellipsoid (default WGS84)
    deg : bool, optional
         degrees input/output  (False: radians in/out)

    Returns
    -------
    geodetic_lat : float or numpy.ndarray of float
         geodetic latiude

    Notes
    -----
    Equations from J. P. Snyder, "Map Projections - A Working Manual",
    US Geological Survey Professional Paper 1395, US Government Printing
    Office, Washington, DC, 1987, pp. 13-18.
    """

    if ell is None:
        ell = Ellipsoid()

    if abs(geocentric_lat) > 90:
        raise ValueError("-90 <= latitude <= 90")

    if deg is True:
        geocentric_lat = radians(geocentric_lat)

    geodetic_lat = atan(tan(geocentric_lat) / (1 - (ell.eccentricity) ** 2))

    if deg is True:
        geodetic_lat = degrees(geodetic_lat)

    return geodetic_lat
