#!/usr/bin/env python3

from setuptools import setup 

with open('README.rst') as f:
	long_description = f.read()
	
setup(name='pymap3d',
      version='0.1',
	  description='3-D coordinate conversion utilities',
	  long_description=long_description,
	  author='Michael Hirsch',
	  author_email='hirsch617@gmail.com',
	  url='https://github.com/scienceopen/pymap3d',
	  install_requires=['astropy','numpy'],
      packages=['pymap3d'],
	  )
	  	  
