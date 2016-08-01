#!/usr/bin/env python
"""
 Ported by Michael Hirsch to Python.
Original work by Joaquim Luis (LGPL), Michael Kleder, et al.
"""
from __future__ import division
from numpy import (absolute, sin, cos, tan, arctan2, atleast_1d,
                   radians, degrees, sign, mod, empty, pi, sqrt, tile, nan)


def vreckon(lat1, lon1, rng, azim, ellipsoid=None):
    """
     VRECKON -  Computes points at a specified azimuth and range in an
                ellipsoidal earth

                   - Using the WGS-84 Earth ellipsoid, travel a given
                     distance along
                     a given azimuth starting at a given initial point,
                     and return the
                     endpoint within a few millimeters of accuracy,
                     using Vincenty's algorithm.

     USAGE:
     lat2,lon2 = vreckon(lat1, lon1, rng, azim)
      Transmits ellipsoid definition (either as [a,b] or [a,f]) as fifth
                argument ELLIPSOID

     VARIABLES:
     lat1 = inital latitude (degrees)
     lon1 = initial longitude (degrees)
     rng  = distance (meters)
                   It can be a scalar or a vector. Latter case computes
                    a series of
                   circles (or arc circles, see azim) centered on X,Y
                    (which are scalars)
     azim = intial azimuth (degrees)
                   "azim" is a scalar or vector
     ellipsoid = two-element ellipsoid vector. Either [a b] or [a f]
                   If omitted, defaults to WGS-84
     lat2, lon2 = second point (degrees)
     a21  = reverse azimuth (degrees), at final point facing back toward
            the  intial point

     Original algorithm source:
     T. Vincenty, "Direct and Inverse Solutions of Geodesics on the
Ellipsoid
     with Application of Nested Equations", Survey Review, vol. 23, no.
176,
     April 1975, pp 88-93.
     Available at: http://www.ngs.noaa.gov/PUBS_LIB/inverse.pdf

     Notes:
     (1) The Vincenty reckoning algorithm was transcribed verbatim into
         JavaScript by Chris Veness. It was modified and translated to
Matlab
         by Michael Kleder. Mr. Veness's website is:

http://www.movable-type.co.uk/scripts/latlong-vincenty-direct.html
     (2) Error correcting code, polar error corrections, WGS84 ellipsoid
         parameters, testing, and comments by Michael Kleder.
     (3) By convention, when starting at a pole, the longitude of the
initial
         point (otherwise meaningless) determines the longitude line
along
         which to traverse, and hence the longitude of the final point.
     (4) The convention noted in (3) above creates a discrepancy with
VDIST
         when the the intial or final point is at a pole. In the VDIST
         function, when traversing from a pole, the azimuth is  0 when
         heading away from the south pole and 180 when heading away from
the
         north pole. In contrast, this VRECKON function uses the azimuth
as
        noted in (3) above when traversing away form a pole.
     (5) In testing, where the traversal subtends no more than 178
degrees,
         this function correctly inverts the VDIST function to within
0.2
         millimeters of distance, 5e-10 degrees of forward azimuth,
         and 5e-10 degrees of reverse azimuth. Precision reduces as test
         points approach antipodal because the precision of VDIST is
reduced
         for nearly antipodal points. (A warning is given by VDIST.)
     (6) Tested but no warranty. Use at your own risk.
     (7) Ver 1.0, Michael Kleder, November 2007
     (8) Ver 2.0, Joaquim Luis, September 2008

     Added ellipsoid and vectorized whenever possible.
     Also, lon2 is always converted to the [-180 180] interval
     Joaquim Luis
     $Id: $
    """

    lat1 = atleast_1d(lat1)
    lon1 = atleast_1d(lon1)
    rng = atleast_1d(rng)
    azim = atleast_1d(azim)

    assert absolute(lat1) <= 90, (
        'VRECKON: Input lat. must be between -90 and 90 deg., inclusive.')
    if lat1.size != 1 and rng.size > 1:
        raise ValueError(
            'VRECKON: Variable ranges are only allowed for a single point.')
    if ellipsoid is not None:  # An ellipsoid vector (with a & b OR a & f)
        a = ellipsoid[0]           # b = ellipsoid(2);
        # Second ellipsoid argument contains flattening instead of minor axis
        if ellipsoid[1] < 1:
            f = ellipsoid[1]
            b = a * (1 - f)
        else:               # Second ellipsoid argument contains minor axis
            f = (a - ellipsoid[1]) / a
    else:   # Supply WGS84 earth ellipsoid axis lengths in meters:
        a = 6378137                 # semimajor axis
        b = 6356752.31424518   # WGS84 earth flattening coefficient definition
        f = (a - b) / a

    lat1 = radians(lat1)            # intial latitude in radians
    lon1 = radians(lon1)            # intial longitude in radians

    # correct for errors at exact poles by adjusting 0.6 millimeters:
    kidx = absolute(pi / 2 - absolute(lat1)) < 1e-10
    lat1[kidx] = sign(lat1[kidx]) * (pi / 2 - (1e-10))

    #  Allow for multiple circles starting from the same point
    if lat1.size == 1 and lon1.size == 1 and rng.size > 1:
        lat1 = tile(lat1, rng.shape)
        lon1 = tile(lon1, rng.shape)

    if rng.size == 1:
        rng = tile(rng, azim.shape)

    if azim.size == 1:
        azim = tile(azim, rng.shape)

    assert rng.size == azim.shape[0], (
        'Range must be a scalar or vector with the same shape as azim.')

    alpha1 = radians(azim)  # inital azimuth in radians
    sinAlpha1 = sin(alpha1)
    cosAlpha1 = cos(alpha1)

    tanU1 = (1 - f) * tan(lat1)
    cosU1 = 1 / sqrt(1 + tanU1 ** 2)
    sinU1 = tanU1 * cosU1
    sigma1 = arctan2(tanU1, cosAlpha1)
    sinAlpha = cosU1 * sinAlpha1
    cosSqAlpha = 1 - sinAlpha * sinAlpha
    uSq = cosSqAlpha * (a**2 - b**2) / b**2
    A = 1 + uSq / 16384 * (4096 + uSq * (-768 + uSq * (320 - 175 * uSq)))
    B = uSq / 1024 * (256 + uSq * (-128 + uSq * (74 - 47 * uSq)))
    sigma = rng / (b * A)
    sigmaP = 2 * pi

    if sigma.size == 1:
        sinSigma = nan
        cosSigma = nan
        cos2SigmaM = nan
        while absolute(sigma - sigmaP) > 1e-12:
            cos2SigmaM = cos(2 * sigma1 + sigma)
            sinSigma = sin(sigma)
            cosSigma = cos(sigma)
            deltaSigma = (B * sinSigma *
                          (cos2SigmaM + B / 4 *
                           (cosSigma * (-1 + 2 * cos2SigmaM * cos2SigmaM) -
                            B / 6 * cos2SigmaM *
                            (-3 + 4 * sinSigma * sinSigma) *
                            (-3 + 4 * cos2SigmaM * cos2SigmaM))))
            sigmaP = sigma
            sigma = rng / (b * A) + deltaSigma
    else:
        # This part is not vectorized
        cos2SigmaM = empty(sigma.shape)
        sinSigma = empty(sigma.shape)
        cosSigma = empty(sigma.shape)

        for k in range(sigma.size):
            while (absolute(sigma[k] - sigmaP) > 1e-12).any():
                cos2SigmaM[k] = cos(2 * sigma1[k] + sigma[k])
                sinSigma[k] = sin(sigma[k])
                cosSigma[k] = cos(sigma[k])
                tmp = 2 * cos2SigmaM[k] * cos2SigmaM[k]
                deltaSigma = (B[k] * sinSigma[k] *
                              (cos2SigmaM[k] + B[k] / 4 *
                               (cosSigma[k] * (-1 + tmp) -
                                B[k] / 6 * cos2SigmaM[k] *
                                (-3 + 4 * sinSigma[k] * sinSigma[k]) *
                                  (-3 + 2 * tmp))))
                sigmaP = sigma[k]
                sigma[k] = rng[k] / (b * A[k]) + deltaSigma

    tmp = sinU1 * sinSigma - cosU1 * cosSigma * cosAlpha1
    lat2 = arctan2(sinU1 * cosSigma + cosU1 * sinSigma * cosAlpha1,
                   (1 - f) * sqrt(sinAlpha * sinAlpha + tmp**2))

    lamb = arctan2(sinSigma * sinAlpha1,
                   cosU1 * cosSigma - sinU1 * sinSigma * cosAlpha1)

    C = f / 16 * cosSqAlpha * (4 + f * (4 - 3 * cosSqAlpha))
    L = lamb - (f * (1 - C) *
                sinAlpha * (sigma + C * sinSigma *
                            (cos2SigmaM + C * cosSigma *
                             (-1 + 2 * cos2SigmaM * cos2SigmaM))))

    lon2 = degrees(lon1 + L)

    # Truncates angles into the [-pi pi] range
    # if (lon2 > pi).any():
    #    lon2 = pi*((absolute(lon2)/pi) -
    #       2*ceil(((absolute(lon2)/pi)-1)/2)) * sign(lon2)

    # lon2 = mod(lon2,360); % follow [0,360] convention
    lon2 = (lon2 + 180) % 360 - 180  # no parenthesis on RHS

    a21 = arctan2(sinAlpha, -tmp)
    a21 = 180. + degrees(a21)  # note direction reversal

    return degrees(lat2), lon2, mod(a21, 360.)
