import functools
from pathlib import Path

import matlab.engine


@functools.cache
def matlab_engine():
    cwd = Path(__file__).parent
    eng = matlab.engine.start_matlab("-nojvm")
    eng.addpath(eng.genpath(str(cwd)), nargout=0)
    return eng
