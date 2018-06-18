#!/usr/bin/env python
# Copyright (c) 2014-2018 Michael Hirsch, Ph.D.

from setuptools import setup, find_packages

install_requires = ['python-dateutil', 'pytz', 'numpy']
tests_require = ['pytest', 'coveralls', 'pyproj', 'flake8', 'mypy']
# %%


setup(name='pymap3d',
      packages=find_packages(),
      version='1.7.0',
      description='pure Python coordinate conversions, following convention of several popular Matlab routines.',
      long_description=open('README.rst').read(),
      long_description_content_type="text/markdown",
      author='Michael Hirsch, Ph.D.',
      url='https://github.com/scivision/pymap3d',
      classifiers=[
                  'Development Status :: 5 - Production/Stable',
                  'Environment :: Console',
                  'Intended Audience :: Science/Research',
                  'License :: OSI Approved :: BSD License',
                  'Operating System :: OS Independent',
                  'Programming Language :: Python :: 3.6',
                  'Programming Language :: Python :: 3.7',
                  'Topic :: Scientific/Engineering :: GIS',
      ],
      install_requires=install_requires,
      tests_require=tests_require,
      extras_require={'tests': tests_require,
                      'full': ['astropy']},
      python_requires='>=3.6',
      )
