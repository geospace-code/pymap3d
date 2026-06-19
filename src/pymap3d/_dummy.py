from datetime import datetime

def eci2ecef(x, y, z, time: datetime, force_non_astropy: bool = False) -> tuple:
    raise ImportError("Numpy required for eci2ecef")

def ecef2eci(x, y, z, time: datetime, force_non_astropy: bool = False) -> tuple:
    raise ImportError("Numpy required for ecef2eci")
