"""
this is from PySatel and gives same result to EIGHT decimal places
"""
from pymap3d.ecef import Ellipsoid
from numpy import sqrt, arctan, arctan2, hypot, degrees


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
