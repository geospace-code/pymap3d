#!/usr/bin/env python
"""
vreckon and vdist are iterative algorithms.
How much does PyPy help over Cpython?

Hmm, PyPy is slower than Cpython..

$ pypy3 tests/benchmark_vincenty.py 10000
2.1160879135131836
0.06056046485900879

$ python tests/benchmark_vincenty.py 10000
0.3325080871582031
0.02107095718383789
"""

import argparse
import shutil
import subprocess
import time
from pathlib import Path

import numpy as np
from pymap3d.vincenty import vdist, vreckon

R = Path(__file__).parent

MATLAB = shutil.which("matlab")

ll0 = (42.0, 82.0)


def bench_vreckon(N: int) -> float:

    sr = np.random.random(N)
    az = np.random.random(N)

    tic = time.monotonic()
    _, _ = vreckon(ll0[0], ll0[1], sr, az)

    return time.monotonic() - tic


def bench_vdist(N: int) -> float:
    lat = np.random.random(N)
    lon = np.random.random(N)

    tic = time.monotonic()
    _, _ = vdist(ll0[0], ll0[1], lat, lon)

    return time.monotonic() - tic


if __name__ == "__main__":
    p = argparse.ArgumentParser()
    p.add_argument("N", help="number of iterations", type=int)
    args = p.parse_args()
    N = args.N

    print(f"vreckon: {bench_vreckon(N):.3f}")
    print(f"vdist: {bench_vdist(N):.3f}")

    if MATLAB:
        print(f"matlab path {R}")
        subprocess.check_call(
            f'matlab -batch "helper_vdist({ll0[0]}, {ll0[1]}, {N})"', text=True, timeout=90, cwd=R
        )
