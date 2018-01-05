#!/usr/bin/env python
install_requires = ['python-dateutil','pytz','numpy','astropy']
tests_require = ['nose','coveralls']
# %%
from setuptools import setup,find_packages

setup(name='pymap3d',
      packages=find_packages(),
      version = '1.2.5',
      description='Python coordinate conversions, following convention of several popular Matlab routines.',
      long_description=open('README.rst').read(),
      author = 'Michael Hirsch, Ph.D.',
      url = 'https://github.com/scivision/pymap3d',
      classifiers=[
      'Intended Audience :: Science/Research',
      'Development Status :: 5 - Production/Stable',
      'License :: OSI Approved :: BSD License',
      'Topic :: Scientific/Engineering :: GIS',
      'Programming Language :: Python :: 3',
      'Programming Language :: Python',
      ],
      install_requires=install_requires,
      tests_require=tests_require,
      extras_require={'tests':tests_require},
      python_requires='>=3.5',  # due to @ matrix mult
	  )

