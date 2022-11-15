"""Minimal class for planetary ellipsoids"""

from __future__ import annotations
from math import sqrt
from dataclasses import dataclass, field
import sys
import warnings
from typing import Dict  # for Python < 3.9

if sys.version_info >= (3, 8):
    from typing import TypedDict
else:
    from typing_extensions import TypedDict


class Model(TypedDict):
    """Ellipsoid parameters"""

    name: str
    a: float
    b: float


@dataclass
class Ellipsoid:
    """
    generate reference ellipsoid parameters

    as everywhere else in pymap3d, distance units are METERS

    Ellipsoid sources
    -----------------

    maupertuis, plessis, everest1830, everest1830m, everest1967,
    airy, bessel, clarke1866, clarke1878, clarke1860, helmert, hayford,
    international1924, krassovsky1940, wgs66, australian, international1967,
    grs67, sa1969, wgs72, iers1989, iers2003:

    - https://en.wikipedia.org/wiki/Earth_ellipsoid#Historical_Earth_ellipsoids
    - https://en.wikibooks.org/wiki/PROJ.4#Spheroid

    wgs84: https://en.wikipedia.org/wiki/World_Geodetic_System#WGS84

    wgs84_mean: https://en.wikipedia.org/wiki/Earth_radius#Mean_radii

    grs80: https://en.wikipedia.org/wiki/GRS_80

    io: https://doi.org/10.1006/icar.1998.5987

    pz90.11: https://structure.mil.ru/files/pz-90.pdf

    gsk2011: https://racurs.ru/downloads/documentation/gost_r_32453-2017.pdf

    mars: https://tharsis.gsfc.nasa.gov/geodesy.html

    mercury, venus, moon, jupiter, saturn, uranus, neptune:

    - https://nssdc.gsfc.nasa.gov/planetary/factsheet/index.html

    feel free to suggest additional ellipsoids
    """

    model: str  # short name
    name: str  # name for printing
    semimajor_axis: float
    semiminor_axis: float
    flattening: float
    thirdflattening: float
    eccentricity: float
    models = field(default_factory=Dict[str, Model])

    def __init__(
        self, semimajor_axis: float, semiminor_axis: float, name: str = "", model: str = ""
    ):
        """
        Ellipsoidal model of world

        Parameters
        ----------
        semimajor_axis : float
            semimajor axis in meters
        semiminor_axis : float
            semiminor axis in meters
        name: str, optional
            Human-friendly name for the ellipsoid
        model: str, optional
            Short name for the ellipsoid
        """

        self.flattening = (semimajor_axis - semiminor_axis) / semimajor_axis
        assert self.flattening >= 0, "flattening must be >= 0"
        self.thirdflattening = (semimajor_axis - semiminor_axis) / (semimajor_axis + semiminor_axis)
        self.eccentricity = sqrt(2 * self.flattening - self.flattening**2)

        self.name = name
        self.model = model
        self.semimajor_axis = semimajor_axis
        self.semiminor_axis = semiminor_axis

    models = {
        # Earth ellipsoids
        "maupertuis": {"name": "Maupertuis (1738)", "a": 6397300.0, "b": 6363806.283},
        "plessis": {"name": "Plessis (1817)", "a": 6376523.0, "b": 6355862.9333},
        "everest1830": {"name": "Everest (1830)", "a": 6377299.365, "b": 6356098.359},
        "everest1830m": {
            "name": "Everest 1830 Modified (1967)",
            "a": 6377304.063,
            "b": 6356103.039,
        },
        "everest1967": {
            "name": "Everest 1830 (1967 Definition)",
            "a": 6377298.556,
            "b": 6356097.55,
        },
        "airy": {"name": "Airy (1830)", "a": 6377563.396, "b": 6356256.909},
        "bessel": {"name": "Bessel (1841)", "a": 6377397.155, "b": 6356078.963},
        "clarke1866": {"name": "Clarke (1866)", "a": 6378206.4, "b": 6356583.8},
        "clarke1878": {"name": "Clarke (1878)", "a": 6378190.0, "b": 6356456.0},
        "clarke1860": {"name": "Clarke (1880)", "a": 6378249.145, "b": 6356514.87},
        "helmert": {"name": "Helmert (1906)", "a": 6378200.0, "b": 6356818.17},
        "hayford": {"name": "Hayford (1910)", "a": 6378388.0, "b": 6356911.946},
        "international1924": {"name": "International (1924)", "a": 6378388.0, "b": 6356911.946},
        "krassovsky1940": {"name": "Krassovsky (1940)", "a": 6378245.0, "b": 6356863.019},
        "wgs66": {"name": "WGS66 (1966)", "a": 6378145.0, "b": 6356759.769},
        "australian": {"name": "Australian National (1966)", "a": 6378160.0, "b": 6356774.719},
        "international1967": {
            "name": "New International (1967)",
            "a": 6378157.5,
            "b": 6356772.2,
        },
        "grs67": {"name": "GRS-67 (1967)", "a": 6378160.0, "b": 6356774.516},
        "sa1969": {"name": "South American (1969)", "a": 6378160.0, "b": 6356774.719},
        "wgs72": {"name": "WGS-72 (1972)", "a": 6378135.0, "b": 6356750.52001609},
        "grs80": {"name": "GRS-80 (1979)", "a": 6378137.0, "b": 6356752.31414036},
        "wgs84": {"name": "WGS-84 (1984)", "a": 6378137.0, "b": 6356752.31424518},
        "wgs84_mean": {"name": "WGS-84 (1984) Mean", "a": 6371008.7714, "b": 6371008.7714},
        "iers1989": {"name": "IERS (1989)", "a": 6378136.0, "b": 6356751.302},
        "pz90.11": {"name": "ПЗ-90 (2011)", "a": 6378136.0, "b": 6356751.3618},
        "iers2003": {"name": "IERS (2003)", "a": 6378136.6, "b": 6356751.9},
        "gsk2011": {"name": "ГСК (2011)", "a": 6378136.5, "b": 6356751.758},
        # Other worlds
        "mercury": {"name": "Mercury", "a": 2440500.0, "b": 2438300.0},
        "venus": {"name": "Venus", "a": 6051800.0, "b": 6051800.0},
        "moon": {"name": "Moon", "a": 1738100.0, "b": 1736000.0},
        "mars": {"name": "Mars", "a": 3396900.0, "b": 3376097.80585952},
        "jupyter": {"name": "Jupiter", "a": 71492000.0, "b": 66770054.3475922},
        "io": {"name": "Io", "a": 1829.7, "b": 1815.8},
        "saturn": {"name": "Saturn", "a": 60268000.0, "b": 54364301.5271271},
        "uranus": {"name": "Uranus", "a": 25559000.0, "b": 24973000.0},
        "neptune": {"name": "Neptune", "a": 24764000.0, "b": 24341000.0},
        "pluto": {"name": "Pluto", "a": 1188000.0, "b": 1188000.0},
    }

    @classmethod
    def from_name(cls, name: str) -> Ellipsoid | None:
        """Create an Ellipsoid from a name."""

        if name not in cls.models:
            warnings.warn(f"{name} model not implemented")
            return None

        return cls(
            cls.models[name]["a"], cls.models[name]["b"], name=cls.models[name]["name"], model=name
        )
