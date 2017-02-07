#!/usr/bin/env python
from setuptools import setup

try:
    import conda.cli
    conda.cli.main('install','--file','requirements.txt')
except Exception as e:
    print(e)
    import pip
    pip.main(['install','-r','requirements.txt'])

setup(name='pymap3d',
      packages=['pymap3d'],
      version = '1.0.1',
      description='Python coordinate conversions, following convention of several popular Matlab routines.',
      author = 'Michael Hirsch, Ph.D.',
      url = 'https://github.com/scienceopen/pymap3d',
      classifiers=[
      'Intended Audience :: Science/Research',
      'Development Status :: 5 - Production/Stable',
      'License :: OSI Approved :: BSD License',
      'Topic :: Scientific/Engineering :: GIS',
      'Programming Language :: Python :: 2.7',
      'Programming Language :: Python :: 3.4',
      'Programming Language :: Python :: 3.5',
      'Programming Language :: Python :: 3.6',
      ],
      install_requires=['geopy'],
	  )

