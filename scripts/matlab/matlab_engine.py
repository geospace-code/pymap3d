import functools
from pathlib import Path
from datetime import datetime
import numpy as np

import matlab.engine


def matlab_engine():
    cwd = Path(__file__).parent
    eng = matlab.engine.start_matlab("-nojvm")
    eng.addpath(eng.genpath(str(cwd)), nargout=0)
    return eng


def mat2np(m):
    """
    Matlab matrix to numpy array
    """
    return np.array(m.tomemoryview()).reshape(m.size, order="F")


def pydt2matdt(eng, utc: datetime):
    """
    Python datetime.dateime to Matlab datetime
    """
    return eng.datetime(utc.year, utc.month, utc.day, utc.hour, utc.minute, utc.second)


@functools.cache
def has_matmap3d(eng) -> bool:
    cwd = Path(__file__).parent
    d = cwd.parents[2] / "matmap3d"
    print(f"Looking in {d} for matmap3d")

    if d.is_dir():
        eng.addpath(str(d), nargout=0)
        return True

    return False


@functools.cache
def has_aerospace(eng: matlab.engine.matlabengine.MatlabEngine) -> bool:
    if eng.matlab_toolbox()["aerospace"]:
        return True

    if not has_matmap3d(eng):
        raise EnvironmentError(f"Matlab {eng.version()} does not have Aerospace Toolbox")

    return False


def has_mapping(eng: matlab.engine.matlabengine.MatlabEngine) -> bool:
    if eng.matlab_toolbox()["mapping"]:
        return True

    if not has_matmap3d(eng):
        raise EnvironmentError(f"Matlab {eng.version()} does not have Mapping Toolbox")

    return False
