#!/usr/bin/env python3

"""
https://github.com/geospace-code/pymap3d/issues/103
"""

from datetime import datetime
import sys

from astropy.time import Time
import astropy
import numpy as np

import pymap3d as pm

try:
    from pymap3d.tests.matlab_engine import matlab_ecef2eci, matlab_engine, has_matmap3d
except ImportError:
    pass

print("Python version:", sys.version)
print("AstroPy version:", astropy.__version__)

lat = 33.6  # deg
lon = 134.3  # deg
alt = 0  # m

dt = datetime(2020, 8, 14, 0, 0, 41)

# %% Astropy time for comparison
astropy_time = Time(dt, scale="utc")
print("----------------------------------------")
print("Astropy Time (UTC):", astropy_time.utc)
print("Julian Date (UTC):", astropy_time.utc.jd)
print("Julian Date (TT):", astropy_time.tt.jd)
print("GMST:", astropy_time.sidereal_time("mean", "greenwich"))

np.set_printoptions(precision=3, suppress=True)

# %% 1. Geodetic to ECEF
ecef = pm.geodetic2ecef(lat, lon, alt)
print("\nECEF Coordinates (meters):")
print(np.array(ecef))

# %% AstroPy ECEF to ECI (J2000)
astropy_eci = np.array(pm.ecef2eci(ecef[0], ecef[1], ecef[2], dt))
print("\nAstroPy: ECI Coordinates (meters):")
print(astropy_eci)

numpy_eci = np.array(pm.ecef2eci(ecef[0], ecef[1], ecef[2], dt, force_non_astropy=True))
print("\nNumpy: ECI Coordinates (meters):")
print(numpy_eci)

print("\nAstroPy - Numpy Difference (ECI meters):", astropy_eci - numpy_eci)

# %% Matlab comparison
if matlab_engine in sys.modules:
    eng = matlab_engine()
    eci_matlab_aerospace = matlab_ecef2eci(eng, False, dt, ecef)
    print("\nMatlab Aerospace Toolbox: ECI Coordinates (meters):")
    print(np.array(eci_matlab_aerospace))
    print(
        "\nAstroPy - Matlab Aerospace Toolbox Difference (ECI meters):",
        astropy_eci - eci_matlab_aerospace,
    )
    print(
        "Numpy - Matlab Aerospace Toolbox Difference (ECI meters):",
        numpy_eci - eci_matlab_aerospace,
    )

    if has_matmap3d(eng):
        eci_matmap3d = matlab_ecef2eci(eng, True, dt, ecef)
        print("\nMatlab Matmap3D: ECI Coordinates (meters):")
        print(np.array(eci_matmap3d))
        print(
            "\nAstroPy - Matlab Matmap3D Difference (ECI meters):",
            astropy_eci - eci_matmap3d,
        )
        print("Numpy - Matlab Matmap3D Difference (ECI meters):", numpy_eci - eci_matmap3d)
