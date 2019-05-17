

class Ellipsoid:
    """
    generate reference ellipsoid parameters

    https://en.wikibooks.org/wiki/PROJ.4#Spheroid

    https://nssdc.gsfc.nasa.gov/planetary/factsheet/index.html

    as everywhere else in this program, distance units are METERS
    """

    # def __init__(self, model: str = 'wgs84'):
    def __init__(self, model='wgs84'):
        """
        feel free to suggest additional ellipsoids

        Parameters
        ----------
        model : str
                name of ellipsoid
        """
        if model == 'wgs84':
            """https://en.wikipedia.org/wiki/World_Geodetic_System#WGS84"""
            self.semimajor_axis = 6378137.
            self.flattening = 1 / 298.2572235630
            self.semiminor_axis = self.semimajor_axis * (1 - self.flattening)
        elif model == 'wgs72':
            self.semimajor_axis = 6378135.
            self.flattening = 1 / 298.26
            self.semiminor_axis = self.semimajor_axis * (1 - self.flattening)
        elif model == 'grs80':
            """https://en.wikipedia.org/wiki/GRS_80"""
            self.semimajor_axis = 6378137.
            self.flattening = 1 / 298.257222100882711243
            self.semiminor_axis = self.semimajor_axis * (1 - self.flattening)
        elif model == 'clarke1866':
            self.semimajor_axis = 6378206.4
            self.semiminor_axis = 6356583.8
            self.flattening = -(self.semiminor_axis / self.semimajor_axis - 1)
        elif model == 'mars':
            """
            https://tharsis.gsfc.nasa.gov/geodesy.html
            """
            self.semimajor_axis = 3396900.
            self.semiminor_axis = 3376097.80585952
            self.flattening = 1 / 163.295274386012
        elif model == 'moon':
            self.semimajor_axis = 1738000.
            self.semiminor_axis = self.semimajor_axis
            self.flattening = 0.
        elif model == 'venus':
            self.semimajor_axis = 6051000.
            self.semiminor_axis = self.semimajor_axis
            self.flattening = 0.
        elif model == 'jupiter':
            self.semimajor_axis = 71492000.
            self.flattening = 1 / 15.415446277169725
            self.flattening = -(self.semiminor_axis / self.semimajor_axis - 1)
        elif model == 'io':
            """
            https://doi.org/10.1006/icar.1998.5987
            """
            self.semimajor_axis = 1829.7
            self.flattening = 1 / 131.633093525179
            self.semiminor_axis = self.semimajor_axis * (1 - self.flattening)
        elif model == 'pluto':
            self.semimajor_axis = 1187000.
            self.semiminor_axis = self.semimajor_axis
            self.flattening = 0.
        else:
            raise NotImplementedError('{} model not implemented, let us know and we will add it (or make a pull request)'.format(model))
