""" isometric latitude, meridian distance """
try:
    from numpy import radians, degrees, cos, arctan2 as atan2, tan, pi, vectorize
except ImportError:
    from math import radians, degrees, cos, atan2, tan, pi

    vectorize = None
import typing
from .ellipsoid import Ellipsoid
from .rcurve import rcurve_parallel
from .rsphere import rsphere_rectifying
from .latitude import geodetic2rectifying, rectifying2geodetic, geodetic2isometric, geodetic2authalic, authalic2geodetic
from .utils import sph2cart, cart2sph

try:
    from numpy.typing import ArrayLike
except ImportError:
    ArrayLike = typing.Any

__all__ = ["loxodrome_inverse", "loxodrome_direct", "meridian_arc", "meridian_dist", "departure", "meanm"]


def meridian_dist(lat: ArrayLike, ell: Ellipsoid = None, deg: bool = True) -> ArrayLike:
    """
    Computes the ground distance on an ellipsoid from the equator to the input latitude.

    Parameters
    ----------
    lat : ArrayLike
        geodetic latitude
    ell : Ellipsoid, optional
         reference ellipsoid (default WGS84)
    deg : bool, optional
         degrees input/output  (False: radians in/out)

    Results
    -------
    dist : ArrayLike
         distance (meters)
    """
    return meridian_arc(0, lat, ell, deg)


def meridian_arc(lat1: ArrayLike, lat2: ArrayLike, ell: Ellipsoid = None, deg: bool = True) -> ArrayLike:
    """
    Computes the ground distance on an ellipsoid between two latitudes.

    Parameters
    ----------
    lat1, lat2 : ArrayLike
        geodetic latitudes
    ell : Ellipsoid, optional
         reference ellipsoid (default WGS84)
    deg : bool, optional
         degrees input/output  (False: radians in/out)

    Results
    -------
    dist : ArrayLike
         distance (meters)
    """

    if deg:
        lat1, lat2 = radians(lat1), radians(lat2)

    rlat1 = geodetic2rectifying(lat1, ell, deg=False)
    rlat2 = geodetic2rectifying(lat2, ell, deg=False)

    return rsphere_rectifying(ell) * abs(rlat2 - rlat1)


def loxodrome_inverse(
    lat1: ArrayLike, lon1: ArrayLike, lat2: ArrayLike, lon2: ArrayLike, ell: Ellipsoid = None, deg: bool = True
) -> typing.Tuple[ArrayLike, ArrayLike]:
    """
    computes the arc length and azimuth of the loxodrome
    between two points on the surface of the reference ellipsoid

    like Matlab distance('rh',...) and azimuth('rh',...)

    Parameters
    ----------

    lat1 : ArrayLike
        geodetic latitude of first point
    lon1 : ArrayLike
        geodetic longitude of first point
    lat2 : ArrayLike
        geodetic latitude of second point
    lon2 : ArrayLike
        geodetic longitude of second point
    ell : Ellipsoid, optional
         reference ellipsoid (default WGS84)
    deg : bool, optional
         degrees input/output  (False: radians in/out)

    Results
    -------

    lox_s : ArrayLike
        distance along loxodrome
    az12 : ArrayLike
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
    if vectorize is not None:
        fun = vectorize(loxodrome_inverse_point)
        lox, az = fun(lat1, lon1, lat2, lon2, ell, deg)
        return lox[()], az[()]
    else:
        return loxodrome_inverse_point(lat1, lon1, lat2, lon2, ell, deg)


def loxodrome_inverse_point(
    lat1: float, lon1: float, lat2: float, lon2: float, ell: Ellipsoid = None, deg: bool = True
) -> typing.Tuple[float, float]:
    if deg:
        lat1, lon1, lat2, lon2 = radians(lat1), radians(lon1), radians(lat2), radians(lon2)

    # compute changes in isometric latitude and longitude between points
    disolat = geodetic2isometric(lat2, deg=False, ell=ell) - geodetic2isometric(lat1, deg=False, ell=ell)
    dlon = lon2 - lon1

    # compute azimuth
    az12 = atan2(dlon, disolat)
    cosaz12 = cos(az12)

    # compute distance along loxodromic curve
    if abs(cosaz12) < 1e-9:  # straight east or west
        dist = departure(lon2, lon1, lat1, ell, deg=False)
    else:
        dist = meridian_arc(lat2, lat1, deg=False, ell=ell) / abs(cos(az12))

    if deg:
        az12 = degrees(az12) % 360.0

    return dist, az12


def loxodrome_direct(
    lat1: ArrayLike, lon1: ArrayLike, rng: ArrayLike, a12: ArrayLike, ell: Ellipsoid = None, deg: bool = True
) -> typing.Tuple[ArrayLike, ArrayLike]:
    """
    Given starting lat, lon with arclength and azimuth, compute final lat, lon

    like Matlab reckon('rh', ...)

    Parameters
    ----------
    lat1 : ArrayLike
        inital geodetic latitude (degrees)
    lon1 : ArrayLike
        initial geodetic longitude (degrees)
    rng : ArrayLike
        ground distance (meters)
    a12 : ArrayLike
        azimuth (degrees) clockwide from north.
    ell : Ellipsoid, optional
          reference ellipsoid
    deg : bool, optional
         degrees input/output  (False: radians in/out)

    Results
    -------
    lat2 : ArrayLike
        final geodetic latitude (degrees)
    lon2 : ArrayLike
        final geodetic longitude (degrees)
    """
    if vectorize is not None:
        fun = vectorize(loxodrome_direct_point)
        lat2, lon2 = fun(lat1, lon1, rng, a12, ell, deg)
        return lat2[()], lon2[()]
    else:
        return loxodrome_direct_point(lat1, lon1, rng, a12, ell, deg)


def loxodrome_direct_point(
    lat1: float, lon1: float, rng: float, a12: float, ell: Ellipsoid = None, deg: bool = True
) -> typing.Tuple[float, float]:

    if deg:
        lat1, lon1, a12 = radians(lat1), radians(lon1), radians(a12)

    if abs(lat1) > pi / 2:
        raise ValueError("-90 <= latitude <= 90")
    if rng < 0:
        raise ValueError("ground distance must be >= 0")

    #   compute rectifying sphere latitude and radius
    reclat = geodetic2rectifying(lat1, ell, deg=False)

    # compute the new points
    cosaz = cos(a12)
    lat2 = reclat + (rng / rsphere_rectifying(ell)) * cosaz  # compute rectifying latitude
    lat2 = rectifying2geodetic(lat2, ell, deg=False)  # transform to geodetic latitude

    newiso = geodetic2isometric(lat2, ell, deg=False)
    iso = geodetic2isometric(lat1, ell, deg=False)

    dlon = tan(a12) * (newiso - iso)
    lon2 = lon1 + dlon

    if deg:
        lat2, lon2 = degrees(lat2), degrees(lon2)

    return lat2, lon2


def departure(lon1: ArrayLike, lon2: ArrayLike, lat: ArrayLike, ell: Ellipsoid = None, deg: bool = True) -> ArrayLike:
    """
    Computes the distance along a specific parallel between two meridians.

    like Matlab departure()

    Parameters
    ----------
    lon1, lon2 : ArrayLike
        geodetic longitudes (degrees)
    lat : ArrayLike
        geodetic latitude (degrees)
    ell : Ellipsoid, optional
          reference ellipsoid
    deg : bool, optional
         degrees input/output  (False: radians in/out)

    Returns
    -------
    dist: ArrayLike
        ground distance (meters)
    """
    if deg:
        lon1, lon2, lat = radians(lon1), radians(lon2), radians(lat)

    return rcurve_parallel(lat, ell, deg=False) * ((lon2 - lon1) % pi)


def meanm(lat: ArrayLike, lon: ArrayLike, ell: Ellipsoid = None, deg: bool = True) -> typing.Tuple[float, float]:
    """
    Computes geographic mean for geographic points on an ellipsoid

    like Matlab meanm()

    Parameters
    ----------
    lat : sequence of float
        geodetic latitude (degrees)
    lon : sequence of float
        geodetic longitude (degrees)
    ell : Ellipsoid, optional
          reference ellipsoid
    deg : bool, optional
         degrees input/output  (False: radians in/out)

    Returns
    -------
    latbar, lonbar: float
        geographic mean latitude, longitude
    """
    if deg:
        lat, lon = radians(lat), radians(lon)

    lat = geodetic2authalic(lat, ell, deg=False)
    x, y, z = sph2cart(lon, lat, 1)
    lonbar, latbar, _ = cart2sph(x.sum(), y.sum(), z.sum())
    latbar = authalic2geodetic(latbar, ell, deg=False)

    if deg:
        latbar, lonbar = degrees(latbar), degrees(lonbar)
    return latbar, lonbar
