import functools
from pathlib import Path
from datetime import datetime

import matlab.engine


@functools.cache
def matlab_engine():
    """
    only cached because used by Pytest in multiple tests
    """
    cwd = Path(__file__).parent
    eng = matlab.engine.start_matlab("-nojvm")
    eng.addpath(eng.genpath(str(cwd)), nargout=0)
    return eng


def pydt2matdt(eng, utc: datetime):
    """
    Python datetime.dateime to Matlab datetime
    """
    return eng.datetime(utc.year, utc.month, utc.day, utc.hour, utc.minute, utc.second)


@functools.cache
def has_matmap3d(eng) -> bool:
    cwd = Path(__file__).parent
    d = cwd.parents[3] / "matmap3d"
    print(f"Looking in {d} for matmap3d")

    if d.is_dir():
        eng.addpath(str(d), nargout=0)
        return True

    return False


@functools.cache
def has_aerospace(eng) -> bool:
    return eng.matlab_toolbox()["aerospace"]


@functools.cache
def has_mapping(eng) -> bool:
    return eng.matlab_toolbox()["mapping"]
