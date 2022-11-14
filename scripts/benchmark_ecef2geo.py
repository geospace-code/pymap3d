#!/usr/bin/env python3
"""
benchmark ecef2geodetic
"""
import argparse
import time

import numpy as np
from pymap3d.ecef import ecef2geodetic

ll0 = (42.0, 82.0)


def bench(N: int) -> float:

    x = np.random.random(N)
    y = np.random.random(N)
    z = np.random.random(N)

    tic = time.monotonic()
    _, _, _ = ecef2geodetic(x, y, z)

    return time.monotonic() - tic


if __name__ == "__main__":
    p = argparse.ArgumentParser()
    p.add_argument("N", type=int)
    args = p.parse_args()
    N = args.N

    print(f"ecef2geodetic: {bench(N):.3f} seconds")
