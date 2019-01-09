#!/usr/bin/env python
"""
Generate test input points for use in the various programs and functions
"""
import numpy as np
from pathlib import Path


# %% AER
az = np.arange(0., 360. + 15, 15.)
el = np.arange(0., 90. + 15, 15.)
rng = np.arange(0., 1. + 0.5, 0.5)

with Path('tests/aer.txt').open('w') as f:
    for a in az:
        for e in el:
            for r in rng:
                f.write(f'{a:5.1f}    {e:5.1f}    {r:5.1f}\n')
