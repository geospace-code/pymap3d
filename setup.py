#!/usr/bin/env python
from setuptools import setup
import subprocess

try:
    import conda.cli
    conda.cli.main('install','--file','requirements.txt')
except Exception as e:
    print(e)

setup(name='pymap3d',
      packages=['pymap3d'],
      install_requires=['geopy'],
	  )

