.. image:: https://zenodo.org/badge/DOI/10.5281/zenodo.213676.svg
   :target: https://doi.org/10.5281/zenodo.213676

.. image:: http://img.shields.io/badge/powered%20by-AstroPy-orange.svg?style=flat
    :target: http://www.astropy.org/

.. image:: https://travis-ci.org/scivision/pymap3d.svg?branch=master
    :target: https://travis-ci.org/scivision/pymap3d
    
.. image:: https://coveralls.io/repos/github/scivision/pymap3d/badge.svg?branch=master
    :target: https://coveralls.io/github/scivision/pymap3d?branch=master

==================================
Python 3-D coordinate conversions
==================================

Python coordinate conversions, following convention of several popular Matlab routines.

Usage
=====
a few quick examples::

   from pymap3d import *

   lat,lon,alt = eci2geodetic(eci,t)
   az,el,range = geodetic2aer(lat,lon,alt,42,-82,0)

Functions
==========
Popular mapping toolbox functions ported to Python include::

  aer2ecef  aer2enu  aer2geodetic  aer2ned
  ecef2aer  ecef2enu  ecef2enuv  ecef2geodetic  ecef2ned  ecef2nedv  ecef2eci
  eci2ecef
  enu2aer  enu2ecef  enu2ecefv  enu2geodetic
  geodetic2aer  geodetic2ecef  geodetic2enu  geodetic2ned
  ned2aer  ned2ecef  ned2ecefv  ned2geodetic
  vreckon vdist
  azel2radec radec2azel


Installation
============
::

  pip install pymap3d


Caveats
=======
Atmospheric effects neglected in all functions not invoking AstroPy.
Planetary perturbations and nutation etc. not fully considered.
