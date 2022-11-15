from pathlib import Path

import matlab.engine


def matlab_aerospace(eng: matlab.engine.matlabengine.MatlabEngine) -> bool:

    if eng.has_aerospace_toolbox():
        has_aero = True
    else:
        has_aero = False

        cwd = Path(__file__).parent
        d = cwd.parents[2] / "matmap3d"
        if d.is_dir():
            eng.addpath(str(d), nargout=0)
        else:
            raise EnvironmentError(f"Matlab {eng.version()} does not have Aerospace Toolbox")

    return has_aero
