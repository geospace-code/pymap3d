#!/usr/bin/env python
from setuptools import setup

setup(name='pymap3d',
      packages=['pymap3d'],
      install_requires=['python-dateutil','pytz','nose','numpy','astropy',
                'geopy'],
	  )

