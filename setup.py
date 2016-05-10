#!/usr/bin/env python3

from setuptools import setup
import subprocess

try:
    subprocess.run(['conda','install','--yes','--file','requirements.txt'])
except Exception as e:
    print('you will need to install packages in requirements.txt  {}'.format(e))


with open('README.rst','r') as f:
	long_description = f.read()

setup(name='pymap3d',
      version='0.1',
	  description='3-D coordinate conversion utilities',
	  long_description=long_description,
	  author='Michael Hirsch',
	  install_requires=['geopy'],
	  url='https://github.com/scienceopen/pymap3d',
          packages=['pymap3d'],
	  )

