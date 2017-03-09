#!/usr/bin/env python
from setuptools import setup

req = ['python-dateutil','pytz','nose','numpy','astropy']

setup(name='pymap3d',
      packages=['pymap3d'],
      version = '1.1.1',
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
      install_requires=req,
	  )

