"""
Transformation of 3D coordinates between geocentric geodetic (latitude,
longitude, height) and geocentric spherical (spherical latitude, longitude,
radius).
"""
from __future__ import annotations

from .ellipsoid import Ellipsoid
from .mathfun import asin, atan2, cbrt, degrees, hypot, power, radians, sin, sqrt
from .utils import sanitize

__all__ = [
    "geodetic2spherical",
    "spherical2geodetic",
]


def geodetic2spherical(
    lat,
    lon,
    alt,
    ell: Ellipsoid = None,
    deg: bool = True,
) -> tuple:
    """
    point transformation from Geodetic of specified ellipsoid (default WGS-84)
    to geocentric spherical of the same ellipsoid

    Parameters
    ----------

    lat
           target geodetic latitude
    lon
           target geodetic longitude
    h
         target altitude above geodetic ellipsoid (meters)
    ell : Ellipsoid, optional
          reference ellipsoid
    deg : bool, optional
          degrees input/output  (False: radians in/out)


    Returns
    -------

    Geocentric spherical (spherical latitude, longitude, radius

    lat
           target spherical latitude
    lon
           target longitude
    radius
         target distance to the geocenter (meters)

    based on:
    Vermeille, H., 2002. Direct transformation from geocentric coordinates to
    geodetic coordinates. Journal of Geodesy. 76. 451-454.
    doi:10.1007/s00190-002-0273-6
    """
    lat, ell = sanitize(lat, ell, deg)
    if deg:
        lon = radians(lon)

    # Pre-compute to avoid repeated trigonometric functions
    sinlat = sin(lat)
    coslat = sqrt(1 - sinlat**2)

    # radius of curvature of the prime vertical section
    N = ell.semimajor_axis**2 / hypot(
        ell.semimajor_axis * coslat,
        ell.semiminor_axis * sinlat,
    )

    # Instead of computing X and Y, we only compute the projection on the XY
    # plane: xy_projection = sqrt( X**2 + Y**2 )
    xy_projection = (alt + N) * coslat
    z_cartesian = (alt + (1 - ell.eccentricity**2) * N) * sinlat
    radius = hypot(xy_projection, z_cartesian)
    slat = asin(z_cartesian / radius)

    if deg:
        slat = degrees(slat)
        lon = degrees(lon)

    return slat, lon, radius


def spherical2geodetic(
    lat,
    lon,
    radius,
    ell: Ellipsoid = None,
    deg: bool = True,
) -> tuple:
    """
    point transformation from geocentric spherical of specified ellipsoid
    (default WGS-84) to geodetic of the same ellipsoid

    Parameters
    ----------
    lat
           target spherical latitude
    lon
           target longitude
    radius
         target distance to the geocenter (meters)
    ell : Ellipsoid, optional
          reference ellipsoid
    deg : bool, optional
          degrees input/output  (False: radians in/out)

    Returns
    -------
    lat
           target geodetic latitude
    lon
           target geodetic longitude
    alt
         target altitude above geodetic ellipsoid (meters)

    based on:
    Vermeille, H., 2002. Direct transformation from geocentric coordinates to
    geodetic coordinates. Journal of Geodesy. 76. 451-454.
    doi:10.1007/s00190-002-0273-6
    """
    lat, ell = sanitize(lat, ell, deg)
    if deg:
        lon = radians(lon)

    # Pre-compute to avoid repeated trigonometric functions
    sinlat = sin(lat)
    coslat = sqrt(1 - sinlat**2)

    Z = radius * sinlat
    p_0 = power(radius, 2) * coslat**2 / ell.semimajor_axis**2
    q_0 = (1 - ell.eccentricity**2) / ell.semimajor_axis**2 * Z**2
    r_0 = (p_0 + q_0 - ell.eccentricity**4) / 6
    s_0 = ell.eccentricity**4 * p_0 * q_0 / 4 / r_0**3
    t_0 = cbrt(1 + s_0 + sqrt(2 * s_0 + s_0**2))
    u_0 = r_0 * (1 + t_0 + 1 / t_0)
    v_0 = sqrt(u_0**2 + q_0 * ell.eccentricity**4)
    w_0 = ell.eccentricity**2 * (u_0 + v_0 - q_0) / 2 / v_0
    k = sqrt(u_0 + v_0 + w_0**2) - w_0
    D = k * radius * coslat / (k + ell.eccentricity**2)
    hypotDZ = hypot(D, Z)

    glat = 2 * atan2(Z, (D + hypotDZ))
    alt = (k + ell.eccentricity**2 - 1) / k * hypotDZ

    if deg:
        glat = degrees(glat)
        lon = degrees(lon)

    return glat, lon, alt
