"""
this is from PySatel and gives same result to EIGHT decimal places
"""
from pymap3d.ecef import Ellipsoid
from numpy import sqrt, arctan, arctan2, hypot, degrees
from typing import Tuple


def cbrt(x):
    if x >= 0:
        return pow(x, 1.0 / 3.0)
    else:
        return -pow(abs(x), 1.0 / 3.0)


def ecef2geodetic(x, y, z, ell=Ellipsoid(), deg=True):
    a = ell.a
    b = ell.b
    esq = 6.69437999014 * 0.001
    e1sq = 6.73949674228 * 0.001
    r = hypot(x, y)
    Esq = a**2 - b**2
    F = 54 * b**2 * z**2
    G = r**2 + (1 - esq) * z**2 - esq * Esq
    C = (esq**2 * F * r**2) / (pow(G, 3))
    S = cbrt(1 + C + sqrt(C**2 + 2 * C))
    P = F / (3 * pow((S + 1 / S + 1), 2) * G**2)
    Q = sqrt(1 + 2 * esq**2 * P)
    r_0 = -(P * esq * r) / (1 + Q) + sqrt(0.5 * a**2 * (1 + 1.0 / Q) -
                                          P * (1 - esq) * z**2 / (Q * (1 + Q)) - 0.5 * P * r**2)
    U = sqrt(pow((r - esq * r_0), 2) + z**2)
    V = sqrt(pow((r - esq * r_0), 2) + (1 - esq) * z**2)
    Z_0 = b**2 * z / (a * V)
    alt = U * (1 - b**2 / (a * V))
    lat = arctan((z + e1sq * Z_0) / r)
    lon = arctan2(y, x)

    if deg:
        return degrees(lat), degrees(lon), alt
    else:
        return lat, lon, alt  # radians


def ecef2geodetic_old(x: float, y: float, z: float,
                      ell: Ellipsoid = None, deg: bool = True) -> Tuple[float, float, float]:
    """
    convert ECEF (meters) to geodetic coordinates

    input
    -----
    x,y,z  [meters] target ECEF location                             [0,Infinity)
    ell    reference ellipsoid
    deg    degrees input/output  (False: radians in/out)

    output
    ------
    lat,lon   (degrees/radians)
    alt  (meters)

    Algorithm is based on
    http://www.astro.uni.torun.pl/~kb/Papers/geod/Geod-BG.htm
    This algorithm provides a converging solution to the latitude equation
    in terms of the parametric or reduced latitude form (v)
    This algorithm provides a uniform solution over all latitudes as it does
    not involve division by cos(phi) or sin(phi)
    """
    if ell is None:
        ell = Ellipsoid()

    ea = ell.a
    eb = ell.b
    rad = hypot(x, y)
# Constant required for Latitude equation
    rho = arctan2(eb * z, ea * rad)
# Constant required for latitude equation
    c = (ea**2 - eb**2) / hypot(ea * rad, eb * z)
# Starter for the Newtons Iteration Method
    vnew = arctan2(ea * z, eb * rad)
# Initializing the parametric latitude
    v = 0
    for _ in range(5):
        v = deepcopy(vnew)
# %% Newtons Method for computing iterations
        vnew = v - ((2 * sin(v - rho) - c * sin(2 * v)) /
                    (2 * (cos(v - rho) - c * cos(2 * v))))

        if allclose(v, vnew):
            break
# %% Computing latitude from the root of the latitude equation
    lat = arctan2(ea * tan(vnew), eb)
    # by inspection
    lon = arctan2(y, x)

    alt = (((rad - ea * cos(vnew)) * cos(lat)) +
           ((z - eb * sin(vnew)) * sin(lat)))

    with np.errstate(invalid='ignore'):
        # NOTE: need np.any() to handle scalar and array cases
        if np.any((lat < -pi / 2) | (lat > pi / 2)):
            raise ValueError('-90 <= lat <= 90')

        if np.any((lon < -pi) | (lon > 2 * pi)):
            raise ValueError('-180 <= lat <= 360')

    if deg:
        return degrees(lat), degrees(lon), alt
    else:
        return lat, lon, alt  # radians
