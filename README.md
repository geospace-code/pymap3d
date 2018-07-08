[![image](https://zenodo.org/badge/DOI/10.5281/zenodo.213676.svg)](https://doi.org/10.5281/zenodo.213676)
[![image](http://joss.theoj.org/papers/10.21105/joss.00580/status.svg)](https://doi.org/10.21105/joss.00580)
[![image](http://img.shields.io/badge/powered%20by-AstroPy-orange.svg?style=flat)](http://www.astropy.org/)
[![image](https://travis-ci.org/scivision/pymap3d.svg?branch=master)](https://travis-ci.org/scivision/pymap3d)
[![image](https://coveralls.io/repos/github/scivision/pymap3d/badge.svg?branch=master)](https://coveralls.io/github/scivision/pymap3d?branch=master)
[![image](https://ci.appveyor.com/api/projects/status/af479t19j66t8x5n?svg=true)](https://ci.appveyor.com/project/scivision/pymap3d)
[![image](https://img.shields.io/pypi/pyversions/pymap3d.svg)](https://pypi.python.org/pypi/pymap3d)
[![image](https://img.shields.io/pypi/format/pymap3d.svg)](https://pypi.python.org/pypi/pymap3d)
[![Maintainability](https://api.codeclimate.com/v1/badges/b6e4b90175e6dbf1b375/maintainability)](https://codeclimate.com/github/scivision/pymap3d/maintainability)
[![PyPi Download stats](http://pepy.tech/badge/pymap3d)](http://pepy.tech/project/pymap3d)

# Python / Matlab / Fortran 3-D coordinate conversions

3-D geographic coordinate conversions, with API similar to popular $1000 Matlab Mapping Toolbox routines for:

-   Python
-   Matlab, GNU Octave
-   modern Fortran 2008 standard ("elemental" functions and subroutines
    for massively parallel computation)

PyMap3D is intended for non-interactive use on massively parallel (HPC)
and embedded systems. 
Includes some relevant 
[Vallado's algorithms](http://www.smad.com/vallado/fortran/fortran.html).

[API docs](https://www.scivision.co/pymap3d)

Why not [PyProj](https://github.com/jswhit/pyproj)?

-   PyMap3D does not require anything beyond pure Python + Numpy.
-   PyMap3D API is similar to Matlab Mapping Toolbox, while PyProj's interface is quite distinct
-   PyMap3D intrinsically handles local coordinate systems such as ENU,
    while for PyProj ENU requires some [additional
    effort](https://github.com/jswhit/pyproj/issues/105).
-   PyProj is oriented towards points on the planet surface, while
    PyMap3D handles points on or above the planet surface equally well,
    particularly important for airborne vehicles and remote sensing.

## Prerequisites

-   Python PyMap3D: Python &gt;= 3.6
    -   [AstroPy](http://www.astropy.org/) (optional): If not present,
        ECI coordinate conversions are not available.
-   Matlab / GNU Octave: under `matlab/`
-   Fortran MapTran: under `fortran/`: Fortran 2018 compliant compiler.  Tested with `gfortran`, `ifort`, PGI and `flang`, where the latter two use Fortran 2003 fallback code.

## Install

The three separate packages are independent, they don't rely on each other.

-   Python PyMap3D:

        pip install pymap3d

    or for the latest development code:

        git clone https://github.com/scivision/pymap3d

        pip install -e .

    One can verify Python functionality after installation by:

        pip install -e .[tests]  
        pytest -v

-   Fortran MapTran:

        cd bin
        
        cmake ..
        
        cmake --build .

    verify Fortran (as well as Python and Matlab/Octave) functionality
    after compiling by:

        ctest -V

-   Matlab/Octave: from within Matlab/Octave:

        addpath('pymap3d/matlab')

    One can verify Matlab code functionality by running:

        tests/Test.m

## Usage

Where consistent with the definition of the functions, all arguments may
be arbitrarily shaped (scalar, N-D array).

### Python

```python
import pymap3d as pm

x,y,z = pm.geodetic2ecef(lat,lon,alt)

az,el,range = pm.geodetic2aer(lat, lon, alt, observer_lat, observer_lon, 0)
```

[Python](https://www.python.org/dev/peps/pep-0448/) 
[argument unpacking](https://docs.python.org/3.6/tutorial/controlflow.html#unpacking-argument-lists)
can be used for compact function arguments with scalars or arbitrarily
shaped N-D arrays:

```python
aer = (az,el,slantrange)
obslla = (obs_lat,obs_lon,obs_alt)

lla = pm.aer2geodetic(*aer,*obslla)
```

where tuple `lla` is comprised of scalar or N-D arrays `(lat,lon,alt)`.

### Matlab / GNU Octave

The syntax is reasonably compatible with the $1000 Matlab Mapping
Toolbox. Under the `matlab/` directory:

```matlab
x,y,z = geodetic2ecef([],lat,lon,alt)

az,el,range = geodetic2aer(lat, lon, alt, observer_lat, observer_lon, observer_alt)
```

### Fortran

The Fortran API under `fortran/` directory is simple like PyMap3D.
Modern Fortran "elemental" procedures throughout enable seamless support
of scalar or array coordinate inputs. Default precision is `real64`, set
at the top of `fortran/maptran.f90`.

```fortran
use maptran

call geodetic2ecef(lat,lon,alt, x,y,z)
call geodetic2aer(lat,lon,alt, observer_lat, observer_lon, observer_alt)
```

### Functions

Popular mapping toolbox functions ported to Python include the
following, where the source coordinate system (before the "2") is
converted to the desired coordinate system:

    aer2ecef  aer2enu  aer2geodetic  aer2ned
    ecef2aer  ecef2enu  ecef2enuv  ecef2geodetic  ecef2ned  ecef2nedv
    ecef2eci  eci2ecef
    enu2aer  enu2ecef   enu2geodetic
    geodetic2aer  geodetic2ecef  geodetic2enu  geodetic2ned
    ned2aer  ned2ecef   ned2geodetic
    azel2radec radec2azel
    vreckon vdist

Abbreviations:

-   [AER: Azimuth, Elevation, Range](https://en.wikipedia.org/wiki/Spherical_coordinate_system)
-   [ECEF: Earth-centered, Earth-fixed](https://en.wikipedia.org/wiki/ECEF)
-   [ECI: Earth-centered Inertial](https://en.wikipedia.org/wiki/Earth-centered_inertial)
-   [ENU: East North Up](https://en.wikipedia.org/wiki/Axes_conventions#Ground_reference_frames:_ENU_and_NED)
-   [NED: North East Down](https://en.wikipedia.org/wiki/North_east_down)
-   [radec: right ascension, declination](https://en.wikipedia.org/wiki/Right_ascension)

### Caveats

-   Atmospheric effects neglected in all functions not invoking AstroPy.
    Would need to update code to add these input parameters (just start
    a GitHub Issue to request).
-   Planetary perturbations and nutation etc. not fully considered.

Matlab / Octave
---------------

The `matlab/` directory contains a subset of the Python conversion
functions, usable from Matlab or GNU Octave. Mathworks currently charges
$1000 for the 
[Matlab Mapping Toolbox](https://www.mathworks.com/products/mapping.html) 
that provides these functions.

-   The full set of Python conversions are accessed from Matlab &ge; R2014b by commands like:

        lla = py.pymap3d.geodetic2ecef(x,y,z)

-   Matlab [documentation](https://www.scivision.co/pymap3d) generated
    by [m2html](https://www.artefact.tk/software/matlab/m2html/).

