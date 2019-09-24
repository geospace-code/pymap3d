"""Minimal class for planetary ellipsoids"""
from math import sqrt


class Ellipsoid:
    """
    generate reference ellipsoid parameters

    https://en.wikibooks.org/wiki/PROJ.4#Spheroid

    https://nssdc.gsfc.nasa.gov/planetary/factsheet/index.html

    as everywhere else in this program, distance units are METERS
    """

    def __init__(self, model: str = "wgs84"):
        """
        feel free to suggest additional ellipsoids

        Parameters
        ----------
        model : str
                name of ellipsoid
        """
        if model == "wgs84":
            """https://en.wikipedia.org/wiki/World_Geodetic_System#WGS84"""
            self.semimajor_axis = 6378137.0
            self.semiminor_axis = 6356752.31424518
        elif model == "wgs72":
            self.semimajor_axis = 6378135.0
            self.semiminor_axis = 6356750.52001609
        elif model == "grs80":
            """https://en.wikipedia.org/wiki/GRS_80"""
            self.semimajor_axis = 6378137.0
            self.semiminor_axis = 6356752.31414036
        elif model == "clarke1866":
            self.semimajor_axis = 6378206.4
            self.semiminor_axis = 6356583.8
        elif model == "mars":
            """
            https://tharsis.gsfc.nasa.gov/geodesy.html
            """
            self.semimajor_axis = 3396900.0
            self.semiminor_axis = 3376097.80585952
        elif model == "moon":
            self.semimajor_axis = 1738000.0
            self.semiminor_axis = self.semimajor_axis
        elif model == "venus":
            self.semimajor_axis = 6051000.0
            self.semiminor_axis = self.semimajor_axis
        elif model == "jupiter":
            self.semimajor_axis = 71492000.0
            self.semiminor_axis = 66770054.3475922
        elif model == "io":
            """
            https://doi.org/10.1006/icar.1998.5987
            """
            self.semimajor_axis = 1829.7
            self.semiminor_axis = 1815.8
        elif model == "pluto":
            self.semimajor_axis = 1187000.0
            self.semiminor_axis = self.semimajor_axis
        else:
            raise NotImplementedError(
                "{} model not implemented, let us know and we will add it (or make a pull request)".format(model)
            )

        self.flattening = (self.semimajor_axis - self.semiminor_axis) / self.semimajor_axis
        self.thirdflattening = (self.semimajor_axis - self.semiminor_axis) / (self.semimajor_axis + self.semiminor_axis)
        self.eccentricity = sqrt(2 * self.flattening - self.flattening ** 2)
