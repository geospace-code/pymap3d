#!/usr/bin/env python
from setuptools import setup
import subprocess

try:
    subprocess.call(['conda','install','--file','requirements.txt'])
except Exception as e:
    pass

setup(name='pymap3d',
      packages=['pymap3d'],
      description='3-D coordinate conversion utilities',
      author='Michael Hirsch',
      install_requires=['pathlib2','geopy'],
      url='https://github.com/scienceopen/pymap3d',
	  )

