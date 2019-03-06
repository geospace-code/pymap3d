# Copyright (c) 2014-2018 Michael Hirsch, Ph.D.
"""
 Ported by Michael Hirsch to Python.
"""
from typing import Tuple
import logging
from numpy import arctan, sqrt, tan, sign, sin, cos, arctan2, arcsin, nan, pi
import numpy as np
from .ecef import Ellipsoid

__all__ = ['vdist', 'vreckon', 'track2']


def vdist(Lat1: float, Lon1: float, Lat2: float, Lon2: float, ell: Ellipsoid = None) -> Tuple[float, float, float]:
    """
    Using the reference ellipsoid, compute the distance between two points
    within a few millimeters of accuracy, compute forward azimuth,
    and compute backward azimuth, all using a vectorized version of
    Vincenty's algorithm:

    ```python
    s,a12,a21 = vdist(lat1, lon1, lat2, lon2, ell)
    ```

    Parameters
    ----------

    Lat1 : float or numpy.ndarray of float
        Geodetic latitude of first point (degrees)
    Lon1 : float or numpy.ndarray of float
        Geodetic longitude of first point (degrees)
    Lat2 : float or numpy.ndarray of float
        Geodetic latitude of second point (degrees)
    Lon2 : float or numpy.ndarray of float
        Geodetic longitude of second point (degrees)
    ell : Ellipsoid, optional
          reference ellipsoid

    Results
    -------

    s : float or numpy.ndarray of float
        distance (meters)
    a12 : float or numpy.ndarray of float
        azimuth (degrees) clockwise from first point to second point (forward)
    a21 : float or numpy.ndarray of float
        azimuth (degrees) clockwise from second point to first point (backward)



    Original algorithm source:
    T. Vincenty, "Direct and Inverse Solutions of Geodesics on the Ellipsoid
    with Application of Nested Equations", Survey Review, vol. 23, no. 176,
    April 1975, pp 88-93.
    Available at: http://www.ngs.noaa.gov/PUBS_LIB/inverse.pdf

    Notes:

     1. lat1,lon1,lat2,lon2 can be any (identical) size/shape. Outputs will have the same size and shape.
     2. Error correcting code, convergence failure traps, antipodal corrections, polar error corrections, WGS84 ellipsoid parameters, testing, and comments: Michael Kleder, 2004.
     3. Azimuth implementation (including quadrant abiguity resolution) and code vectorization, Michael Kleder, Sep 2005.
     4. Vectorization is convergence sensitive; that is, quantities which have already converged to within tolerance are not recomputed during subsequent iterations (while other quantities are still converging).
     5. Vincenty describes his distance algorithm as precise to within 0.01 millimeters, subject to the ellipsoidal model.
     6. For distance calculations, essentially antipodal points are treated as exactly antipodal, potentially reducing accuracy slightly.
     7. Distance failures for points exactly at the poles are eliminated by moving the points by 0.6 millimeters.
     8. The Vincenty distance algorithm was transcribed verbatim by Peter Cederholm, August 12, 2003. It was modified and translated to English by Michael Kleder. Mr. Cederholm's website is http://www.plan.aau.dk/~pce/
     9. Distances agree with the Mapping Toolbox, version 2.2 (R14SP3) with a max relative difference of about 5e-9, except when the two points are nearly antipodal, and except when one point is near the equator and the two longitudes are nearly 180 degrees apart. This function (vdist) is more accurate in such cases.
        For example, note this difference (as of this writing):

        ```python
        vdist(0.2,305,15,125)
        ```

        > 18322827.0131551

        ```python
        distance(0.2,305,15,125,[6378137 0.08181919])
        ```

        > 0

     10. Azimuths FROM the north pole (either forward starting at the north pole or backward when ending at the north pole) are set to 180 degrees by convention.
         Azimuths FROM the south pole are set to 0 degrees by convention.
     11. Azimuths agree with the Mapping Toolbox, version 2.2 (R14SP3) to within about a hundred-thousandth of a degree, except when traversing to or from a pole, where the convention for this function is described in (10), and except in the cases noted above in (9).
     12. No warranties; use at your own risk.
    """
    if ell is None:
        ell = Ellipsoid()
# %% prepare inputs
    lat1 = np.atleast_1d(Lat1)
    lat2 = np.atleast_1d(Lat2)
    lon1 = np.atleast_1d(Lon1)
    lon2 = np.atleast_1d(Lon2)

    assert lat1.shape == lon1.shape and lat2.shape == lon2.shape

    if lat1.shape != lat2.shape:
        if lat1.size == 1:
            lat1 = np.broadcast_to(lat1, lat2.shape)
            lon1 = np.broadcast_to(lon1, lon2.shape)

        if lat2.size == 1:
            lat2 = np.broadcast_to(lat2, lat1.shape)
            lon2 = np.broadcast_to(lon2, lon1.shape)
# %% Input check:
    if ((abs(lat1) > 90) | (abs(lat2) > 90)).any():
        raise ValueError('Input latitudes must be in [-90, 90] degrees.')
# %% Supply WGS84 earth ellipsoid axis lengths in meters:
    a = ell.a
    b = ell.b
# %% preserve true input latitudes:
    lat1tr = lat1.copy()
    lat2tr = lat2.copy()
# %% convert inputs in degrees to np.radians:
    lat1 = np.radians(lat1)
    lon1 = np.radians(lon1)
    lat2 = np.radians(lat2)
    lon2 = np.radians(lon2)
# %% correct for errors at exact poles by adjusting 0.6 millimeters:
    kidx = abs(pi / 2 - abs(lat1)) < 1e-10
    if kidx.any():
        lat1[kidx] = sign(lat1[kidx]) * (pi / 2 - (1e-10))

    kidx = abs(pi / 2 - abs(lat2)) < 1e-10
    if kidx.any():
        lat2[kidx] = sign(lat2[kidx]) * (pi / 2 - (1e-10))

    f = (a - b) / a
    U1 = arctan((1 - f) * tan(lat1))
    U2 = arctan((1 - f) * tan(lat2))
    lon1 = lon1 % (2 * pi)
    lon2 = lon2 % (2 * pi)
    L = abs(lon2 - lon1)
    kidx = L > pi
    if kidx.any():
        L[kidx] = 2 * pi - L[kidx]

    lamb = L.copy()  # NOTE: program will fail without copy!
    lambdaold = np.zeros(lat1.shape)
    itercount = 0
    notdone = np.ones(lat1.shape, dtype=bool)
    alpha = np.zeros(lat1.shape)
    sigma = np.zeros(lat1.shape)
    cos2sigmam = np.zeros(lat1.shape)
    C = np.zeros(lat1.shape)
    warninggiven = False
    sinsigma = np.empty(notdone.shape)
    cossigma = np.empty(notdone.shape)
    while notdone.any():  # force at least one execution
        itercount += 1
        if itercount > 50:
            if not warninggiven:
                logging.warning('Essentially antipodal points--precision may be reduced slightly.')

            lamb[notdone] = pi
            break

        lambdaold[notdone] = lamb[notdone]

        sinsigma[notdone] = sqrt(
            (cos(U2[notdone]) * sin(lamb[notdone]))**2 +
            (cos(U1[notdone]) * sin(U2[notdone]) - sin(U1[notdone]) *
             cos(U2[notdone]) * cos(lamb[notdone]))**2)

        cossigma[notdone] = (sin(U1[notdone]) * sin(U2[notdone]) +
                             cos(U1[notdone]) * cos(U2[notdone]) *
                             cos(lamb[notdone]))
        # eliminate rare imaginary portions at limit of numerical precision:
        sinsigma[notdone] = sinsigma[notdone].real
        cossigma[notdone] = cossigma[notdone].real

        sigma[notdone] = arctan2(sinsigma[notdone], cossigma[notdone])

        alpha[notdone] = (arcsin(cos(U1[notdone]) * cos(U2[notdone]) *
                                 sin(lamb[notdone]) / sin(sigma[notdone])))

        cos2sigmam[notdone] = (cos(sigma[notdone]) - 2 * sin(U1[notdone]) *
                               sin(U2[notdone]) / cos(alpha[notdone])**2)

        C[notdone] = (f / 16 * cos(alpha[notdone])**2 *
                      (4 + f * (4 - 3 * cos(alpha[notdone])**2)))

        lamb[notdone] = (L[notdone] + (1 - C[notdone]) * f * sin(alpha[notdone]) *
                         (sigma[notdone] + C[notdone] * sin(sigma[notdone]) *
                          (cos2sigmam[notdone] + C[notdone] * cos(sigma[notdone]) *
                           (-1 + 2. * cos2sigmam[notdone]**2))))
        # print(f'then, lambda(21752) = {lamb[21752],20})
        # correct for convergence failure for essentially antipodal points
        if (lamb[notdone] > pi).any():
            logging.warning('Essentially antipodal points encountered. Precision may be reduced slightly.')
            warninggiven = True
            lambdaold[lamb > pi] = pi
            lamb[lamb > pi] = pi

        notdone = abs(lamb - lambdaold) > 1e-12

    u2 = cos(alpha)**2 * (a**2 - b**2) / b**2
    A = 1 + u2 / 16384 * (4096 + u2 * (-768 + u2 * (320 - 175 * u2)))
    B = u2 / 1024 * (256 + u2 * (-128 + u2 * (74 - 47 * u2)))
    deltasigma = (B * sin(sigma) *
                  (cos2sigmam + B / 4 * (cos(sigma) * (-1 + 2 * cos2sigmam**2) -
                                         B / 6 * cos2sigmam * (-3 + 4 * sin(sigma)**2) * (-3 + 4 * cos2sigmam**2))))

    dist_m = (b * A * (sigma - deltasigma))

# %% From point #1 to point #2
    # correct sign of lambda for azimuth calcs:
    lamb = abs(lamb)
    kidx = sign(sin(lon2 - lon1)) * sign(sin(lamb)) < 0
    lamb[kidx] = -lamb[kidx]
    numer = cos(U2) * sin(lamb)
    denom = cos(U1) * sin(U2) - sin(U1) * cos(U2) * cos(lamb)
    a12 = arctan2(numer, denom)
    kidx = a12 < 0
    a12[kidx] = a12[kidx] + 2 * pi
    # %% from poles
    a12[lat1tr <= -90] = 0
    a12[lat1tr >= 90] = pi
    az = np.degrees(a12)

# %% From point #2 to point #1
    # correct sign of lambda for azimuth calcs:
    lamb = abs(lamb)
    kidx = sign(sin(lon1 - lon2)) * sign(sin(lamb)) < 0
    lamb[kidx] = -lamb[kidx]
    numer = cos(U1) * sin(lamb)
    denom = sin(U1) * cos(U2) - cos(U1) * sin(U2) * cos(lamb)
    a21 = arctan2(numer, denom)
    kidx = a21 < 0
    a21[kidx] = a21[kidx] + 2 * pi
    # %% backwards from poles:
    a21[lat2tr >= 90] = pi
    a21[lat2tr <= -90] = 0.
    backaz = np.degrees(a21)

    return dist_m.squeeze()[()], az.squeeze()[()], backaz.squeeze()[()]


def vreckon(Lat1: float, Lon1: float, Rng: float, Azim: float,
            ell: Ellipsoid = None) -> Tuple[float, float, float]:
    """
    This is the Vincenty "forward" solution.

    Computes points at a specified azimuth and range in an ellipsoidal earth.
    Using the reference ellipsoid, travel a given distance along a given azimuth starting at a given initial point,
    and return the endpoint within a few millimeters of accuracy, using Vincenty's algorithm.

    ```python
    lat2, lon2 = vreckon(lat1, lon1, rng, azim)
    ```

    Parameters
    ----------

    Lat1 : float
        inital geodetic latitude (degrees)
    Lon1 : float
        initial geodetic longitude (degrees)
    Rng : float
        ground distance (meters)
    Azim : float
        intial azimuth (degrees) clockwide from north.
    ell : Ellipsoid, optional
          reference ellipsoid

    Results
    -------

    Lat2 : float
        final geodetic latitude (degrees)
    Lon2 : float
        final geodetic longitude (degrees)
    a21 : float
        reverse azimuth (degrees), at final point facing back toward the intial point


    Original algorithm: T. Vincenty, "Direct and Inverse Solutions of Geodesics on the Ellipsoid with Application of Nested Equations", Survey Review, vol. 23, no. 176, April 1975, pp 88-93. http://www.ngs.noaa.gov/PUBS_LIB/inverse.pdf

    Notes:

      1. The Vincenty reckoning algorithm was transcribed verbatim into JavaScript by Chris Veness.
         It was modified and translated to Matlab by Michael Kleder.
         Mr. Veness's website is: http://www.movable-type.co.uk/scripts/latlong-vincenty-direct.html
      2. Error correcting code, polar error corrections, WGS84 ellipsoid parameters, testing, and comments by Michael Kleder.
      3. By convention, when starting at a pole, the longitude of the initial point (otherwise meaningless) determines the longitude line along which to traverse, and hence the longitude of the final point.
      4. The convention noted in (3) above creates a discrepancy with VDIST when the the intial or final point is at a pole. In the VDIST
         function, when traversing from a pole, the azimuth is  0 when
         heading away from the south pole and 180 when heading away from the north pole. In contrast, this VRECKON function uses the azimuth as noted in (3) above when traversing away form a pole.
      5. In testing, where the traversal subtends no more than 178 degrees, this function correctly inverts the VDIST function to within 0.2 millimeters of distance, 5e-10 degrees of forward azimuth,
         and 5e-10 degrees of reverse azimuth. Precision reduces as test
         points approach antipodal because the precision of VDIST is reduced
         for nearly antipodal points. (A warning is given by VDIST.)
      6. Tested but no warranty. Use at your own risk.
      7. Ver 1.0, Michael Kleder, November 2007. Ver 2.0, Joaquim Luis, September 2008

    Added ellipsoid and vectorized whenever possible. Also, lon2 is always converted to the [-180 180] interval.
    Joaquim Luis
    """

    lat1 = np.atleast_1d(Lat1)
    lon1 = np.atleast_1d(Lon1)
    rng = np.atleast_1d(Rng)
    azim = np.atleast_1d(Azim)

    if rng.ndim != 1 or azim.ndim != 1:
        raise ValueError('Range and azimuth must be scalar or vector')

    if abs(lat1) > 90:
        raise ValueError('VRECKON: Input lat. must be between -90 and 90 deg., inclusive.')

    if lat1.size > 1 and rng.size > 1:
        raise ValueError('VRECKON: Variable ranges are only allowed for a single point.')

    if ell is not None:
        a = ell.a
        b = ell.b
        f = ell.f
    else:   # Supply WGS84 earth ellipsoid axis lengths in meters:
        a = 6378137                 # semimajor axis
        b = 6356752.31424518   # WGS84 earth flattening coefficient definition
        f = (a - b) / a

    lat1 = np.radians(lat1)            # intial latitude in radians
    lon1 = np.radians(lon1)            # intial longitude in radians

    # correct for errors at exact poles by adjusting 0.6 millimeters:
    kidx = abs(pi / 2 - abs(lat1)) < 1e-10
    lat1[kidx] = sign(lat1[kidx]) * (pi / 2 - (1e-10))

    if rng.size != azim.size and rng.size == 1:
        rng = np.broadcast_to(rng, azim.size)

    alpha1 = np.radians(azim)  # inital azimuth in radians
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
        while abs(sigma - sigmaP) > 1e-12:
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
        cos2SigmaM = np.empty(sigma.size)
        sinSigma = np.empty(sigma.size)
        cosSigma = np.empty(sigma.size)

        for k in range(sigma.size):
            while (abs(sigma[k] - sigmaP) > 1e-12).any():
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

    lon2 = np.degrees(lon1 + L)

    # Truncates angles into the [-pi pi] range
    # if (lon2 > pi).any():
    #    lon2 = pi*((absolute(lon2)/pi) -
    #       2*ceil(((absolute(lon2)/pi)-1)/2)) * sign(lon2)

    # lon2 = mod(lon2,360); % follow [0,360] convention
    lon2 = (lon2 + 180) % 360 - 180  # no parenthesis on RHS

    a21 = arctan2(sinAlpha, -tmp)
    a21 = 180. + np.degrees(a21)  # note direction reversal

    return np.degrees(lat2).squeeze()[()], lon2.squeeze()[()], a21.squeeze()[()] % 360.


def track2(lat1: float, lon1: float, lat2: float, lon2: float,
           ell: Ellipsoid = None, npts: int = 100, deg: bool = True):
    """
    computes great circle tracks starting at the point lat1, lon1 and ending at lat2, lon2

    Parameters
    ----------

    Lat1 : float or numpy.ndarray of float
        Geodetic latitude of first point (degrees)
    Lon1 : float or numpy.ndarray of float
        Geodetic longitude of first point (degrees)
    Lat2 : float or numpy.ndarray of float
        Geodetic latitude of second point (degrees)
    Lon2 : float or numpy.ndarray of float
        Geodetic longitude of second point (degrees)
    ell : Ellipsoid, optional
          reference ellipsoid
    npts : int, optional
        number of points (default is 100)
    deg : bool, optional
        degrees input/output  (False: radians in/out)

    Results
    -------

    lats : numpy.ndarray of float
        latitudes of points along track
    lons : numpy.ndarray of float
        longitudes of points along track

    Based on code posted to the GMT mailing list in Dec 1999 by Jim Levens and by Jeff Whitaker <jeffrey.s.whitaker@noaa.gov>
    """

    if ell is None:
        ell = Ellipsoid()

    if npts <= 1:
        raise ValueError('npts must be greater than 1')

    if npts == 2:
        return [lat1, lat2], [lon1, lon2]

    if deg is True:
        rlat1, rlon1, rlat2, rlon2 = np.radians([lat1, lon1, lat2, lon2])
    else:
        rlat1, rlon1, rlat2, rlon2 = lat1, lon1, lat2, lon2

    gcarclen = 2. * np.arcsin(np.sqrt((np.sin((rlat1 - rlat2) / 2))**2 +
                                      np.cos(rlat1) * np.cos(rlat2) * (np.sin((rlon1 - rlon2) / 2))**2))
    # check to see if points are antipodal (if so, route is undefined).
    if np.allclose(gcarclen, pi):
        raise ValueError('cannot compute intermediate points on a great circle whose endpoints are antipodal')

    distance, azimuth, _ = vdist(lat1, lon1, lat2, lon2)
    incdist = distance / (npts - 1)

    latpt = lat1
    lonpt = lon1
    lons = [lonpt]
    lats = [latpt]
    for n in range(npts - 2):
        latptnew, lonptnew, _ = vreckon(latpt, lonpt, incdist, azimuth)
        _, azimuth, _ = vdist(latptnew, lonptnew, lat2, lon2, ell=ell)
        lats.append(latptnew)
        lons.append(lonptnew)
        latpt = latptnew
        lonpt = lonptnew
    lons.append(lon2)
    lats.append(lat2)

    if not deg:
        lats = np.radians(lats)
        lons = np.radians(lons)

    return lats, lons
