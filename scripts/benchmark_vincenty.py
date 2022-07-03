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
import time
from pymap3d.vincenty import vreckon, vdist
import numpy as np
import argparse
import subprocess
import shutil

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
    p.add_argument("N", type=int)
    p = p.parse_args()
    N = p.N

    print(f"vreckon: {bench_vreckon(N):.3f}")
    print(f"vdist: {bench_vdist(N):.3f}")

    if MATLAB:
        subprocess.check_call(
            f'matlab -batch "f = @() distance({ll0[0]}, {ll0[1]}, rand({N},1), rand({N},1));  t = timeit(f); disp(t)"',
            timeout=45,
        )
