.. image:: https://zenodo.org/badge/DOI/10.5281/zenodo.213676.svg
   :target: https://doi.org/10.5281/zenodo.213676
   
.. image:: http://joss.theoj.org/papers/10.21105/joss.00580/status.svg
    :target: https://doi.org/10.21105/joss.00580

.. image:: http://img.shields.io/badge/powered%20by-AstroPy-orange.svg?style=flat
    :target: http://www.astropy.org/

.. image:: https://travis-ci.org/scivision/pymap3d.svg?branch=master
    :target: https://travis-ci.org/scivision/pymap3d

.. image:: https://coveralls.io/repos/github/scivision/pymap3d/badge.svg?branch=master
    :target: https://coveralls.io/github/scivision/pymap3d?branch=master
    
.. image:: https://ci.appveyor.com/api/projects/status/af479t19j66t8x5n?svg=true
    :target: https://ci.appveyor.com/project/scivision/pymap3d

.. image:: https://img.shields.io/pypi/pyversions/pymap3d.svg
  :target: https://pypi.python.org/pypi/pymap3d
  :alt: Python versions (PyPI)

.. image::  https://img.shields.io/pypi/format/pymap3d.svg
  :target: https://pypi.python.org/pypi/pymap3d
  :alt: Distribution format (PyPI)

.. image:: https://api.codeclimate.com/v1/badges/b6e4b90175e6dbf1b375/maintainability
   :target: https://codeclimate.com/github/scivision/pymap3d/maintainability
   :alt: Maintainability

====================================================
Python / Matlab / Fortran 3-D coordinate conversions
====================================================

3-D geographic coordinate conversions, with API similar to popular $1000 Matlab Mapping Toolbox routines for:

* Python
* Matlab, GNU Octave
* modern Fortran ("elemental" functions and subroutines for massively parallel computation)

PyMap3D is intended for non-interactive use on massively parallel (HPC) and embedded systems.
Includes some relevant `Vallado's algorithms <http://www.smad.com/vallado/fortran/fortran.html>`_.

:API docs: https://www.scivision.co/pymap3d

For those not having:

* `AstroPy <http://www.astropy.org/>`_: lower accuracy fallback functions are automatically used.
* Numpy: without Numpy, scalar inputs are handled with pure Python builtins.

Why not `PyProj <https://github.com/jswhit/pyproj>`_?

* PyMap3D does not require anything beyond pure Python.
* PyMap3D API is virtually identical to Matlab Mapping Toolbox, while PyProj's interface is quite distinct
* PyMap3D intrinsically handles local coordinate systems such as ENU, while for PyProj ENU requires some `additional effort <https://github.com/jswhit/pyproj/issues/105>`_.
* PyProj is oriented towards points on the planet surface, while PyMap3D handles points on or above the planet surface equally well, particularly import for airborne vehicles and remote sensing.

.. contents::


Prerequisites
=============

* Python PyMap3D:  any of Python 2.6, 2.7, 3.4, 3.5, 3.6, but only Python >= 3.5 is tested regularly.
  * Numpy (optional): if not present, function from ``math`` are automatically used, and if Numpy is present, it is automatically used.
  * `AstroPy <http://www.astropy.org/>`_  (optional): If not present, ECI coordinate conversions are not available.
* Matlab / GNU Octave: under ``matlab/``
* Fortran MapTran: under ``fortran/``:  any Fortran compiler (tested with ``gfortran``)

Install
=======
The three separate packages are independent, they don't rely on each other.

* Python PyMap3D::

      pip install pymap3d

  or for the latest development code::

      git clone https://github.com/scivision/pymap3d

      pip install -e .
      
      
  One can verify Python functionality after installation by::
  
      pip install -e .[tests]  
      pytest -v

* Fortran MapTran::

    cd bin
    cmake ..
    make
    
  verify Fortran (as well as Python and Matlab/Octave) functionality after compiling by::
  
    make test

* Matlab/Octave: from within Matlab/Octave::

    addpath(pymap3d/matlab)
    
  One can verify Matlab code functionality by running::
  
      tests/Test.m


Usage
=====

Where consistent with the definition of the functions, all arguments may be arbitrarily shaped (scalar, N-D array).

Python
------

.. code:: python

   import pymap3d as pm

   x,y,z = pm.geodetic2ecef(lat,lon,alt)

   az,el,range = pm.geodetic2aer(lat, lon, alt, observer_lat, observer_lon, 0)
   
`Python >= 3.5 <https://www.python.org/dev/peps/pep-0448/>`_
`argument unpacking <https://docs.python.org/3.6/tutorial/controlflow.html#unpacking-argument-lists>`_ 
can be used for compact function arguments with scalars or arbitrarily shaped N-D arrays:

.. code:: python

    aer = (az,el,slantrange)
    obslla = (obs_lat,obs_lon,obs_alt)
    
    lla = pm.aer2geodetic(*aer,*obslla)
    
where tuple ``lla`` is comprised of scalar or N-D arrays ``(lat,lon,alt)``.



Matlab / GNU Octave
-------------------
The syntax is reasonably compatible with the $1000 Matlab Mapping Toolbox.
Under the ``matlab/`` directory:

.. code:: matlab

   x,y,z = geodetic2ecef([],lat,lon,alt)

   az,el,range = geodetic2aer(lat, lon, alt, observer_lat, observer_lon, observer_alt)


Fortran
-------
The Fortran API under ``fortran/`` directory is simple like PyMap3D.
Modern Fortran "elemental" procedures throughout enable seamless support of scalar or array coordinate inputs.
Default precision is ``real64``, set at the top of ``fortran/maptran.f90``.

.. code:: fortran

    use maptran

    call geodetic2ecef(lat,lon,alt, x,y,z)
    call geodetic2aer(lat,lon,alt, observer_lat, observer_lon, observer_alt)




Functions
---------
Popular mapping toolbox functions ported to Python include the following, where the source coordinate system (before the "2") is converted to the desired coordinate system::

  aer2ecef  aer2enu  aer2geodetic  aer2ned
  ecef2aer  ecef2enu  ecef2enuv  ecef2geodetic  ecef2ned  ecef2nedv
  ecef2eci  eci2ecef
  enu2aer  enu2ecef   enu2geodetic
  geodetic2aer  geodetic2ecef  geodetic2enu  geodetic2ned
  ned2aer  ned2ecef   ned2geodetic
  azel2radec radec2azel
  vreckon vdist

Abbreviations:

* `AER: Azimuth, Elevation, Range <https://en.wikipedia.org/wiki/Spherical_coordinate_system>`_
* `ECEF: Earth-centered, Earth-fixed <https://en.wikipedia.org/wiki/ECEF>`_
* `ECI: Earth-centered Inertial <https://en.wikipedia.org/wiki/Earth-centered_inertial>`_
* `ENU: East North Up <https://en.wikipedia.org/wiki/Axes_conventions#Ground_reference_frames:_ENU_and_NED>`_
* `NED: North East Down <https://en.wikipedia.org/wiki/North_east_down>`_
* `radec: right ascension, declination <https://en.wikipedia.org/wiki/Right_ascension>`_


Caveats
-------

* Atmospheric effects neglected in all functions not invoking AstroPy. Would need to update code to add these input parameters (just start a GitHub Issue to request).
* Planetary perturbations and nutation etc. not fully considered.


Matlab / Octave
===============

The ``matlab/`` directory contains a subset of the Python conversion functions, usable from Matlab or GNU Octave.
Mathworks currently charges $1000 for the `Matlab Mapping Toolbox <https://www.mathworks.com/products/mapping.html>`_ that provides these functions.

* The full set of Python conversions can be accessed from Matlab >= R2014b by commands like::

    lla = py.pymap3d.geodetic2ecef(x,y,z)

* Matlab `documentation <https://www.scivision.co/pymap3d>`_ generated by `m2html <https://www.artefact.tk/software/matlab/m2html/>`_.

