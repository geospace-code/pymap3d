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
from time import time
from pymap3d.vincenty import vreckon, vdist
import numpy as np
from argparse import ArgumentParser

ll0 = (42., 82.)


def bench_vreckon(N: int) -> float:

    sr = np.random.random(N)
    az = np.random.random(N)

    tic = time()
    a, b, c = vreckon(*ll0, sr, az)

    return time() - tic


def bench_vdist(N: int) -> float:
    lat = np.random.random(N)
    lon = np.random.random(N)

    tic = time()
    asr, aaz, aa21 = vdist(*ll0, lat, lon)

    return time() - tic


if __name__ == '__main__':
    p = ArgumentParser()
    p.add_argument('N', type=int)
    p = p.parse_args()

    print(bench_vreckon(p.N))
    print(bench_vdist(p.N))
