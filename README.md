# Python 3-D coordinate conversions

[![image](https://zenodo.org/badge/DOI/10.5281/zenodo.213676.svg)](https://doi.org/10.5281/zenodo.213676)
[![image](http://joss.theoj.org/papers/10.21105/joss.00580/status.svg)](https://doi.org/10.21105/joss.00580)
[![codecov](https://codecov.io/gh/geospace-code/pymap3d/branch/main/graph/badge.svg?token=DFWBW6TKNr)](https://codecov.io/gh/geospace-code/pymap3d)
![Actions Status](https://github.com/geospace-code/pymap3d/workflows/ci/badge.svg)
![Actions Status](https://github.com/geospace-code/pymap3d/workflows/ci_stdlib_only/badge.svg)
[![image](https://img.shields.io/pypi/pyversions/pymap3d.svg)](https://pypi.python.org/pypi/pymap3d)
[![PyPi Download stats](http://pepy.tech/badge/pymap3d)](http://pepy.tech/project/pymap3d)

Pure Python (no prerequistes beyond Python itself) 3-D geographic coordinate conversions and geodesy.
Function syntax is roughly similar to Matlab Mapping Toolbox.
PyMap3D is intended for non-interactive use on massively parallel (HPC) and embedded systems.

[API docs](https://geospace-code.github.io/pymap3d/)

Thanks to our [contributors](./.github/contributors.md).

## Similar toolboxes in other code languages

* [Matlab, GNU Octave](https://github.com/geospace-code/matmap3d)
* [Fortran](https://github.com/geospace-code/maptran3d)
* [Rust](https://github.com/gberrante/map_3d)
* [C++](https://github.com/ClancyWalters/cppmap3d)

## Prerequisites

Numpy and AstroPy are optional.
Algorithms from Vallado and Meeus are used if AstroPy is not present.

## Install

```sh
python3 -m pip install pymap3d
```

or for the latest development code:

```sh
git clone https://github.com/geospace-code/pymap3d

pip install -e pymap3d
```

One can verify Python functionality after installation by:

```sh
pytest pymap3d
```

## Usage

Where consistent with the definition of the functions, all arguments may
be arbitrarily shaped (scalar, N-D array).

```python
import pymap3d as pm

x,y,z = pm.geodetic2ecef(lat,lon,alt)

az,el,range = pm.geodetic2aer(lat, lon, alt, observer_lat, observer_lon, 0)
```

[Python](https://www.python.org/dev/peps/pep-0448/)
[argument unpacking](https://docs.python.org/3/tutorial/controlflow.html#unpacking-argument-lists)
can be used for compact function arguments with scalars or arbitrarily
shaped N-D arrays:

```python
aer = (az,el,slantrange)
obslla = (obs_lat ,obs_lon, obs_alt)

lla = pm.aer2geodetic(*aer, *obslla)
```

where tuple `lla` is comprised of scalar or N-D arrays `(lat,lon,alt)`.

Example scripts are in the [examples](./Examples) directory.

Native Python float is typically [64 bit](https://docs.python.org/3/library/stdtypes.html#typesnumeric).
Numpy can select real precision bits: 32, 64, 128, etc.

There is also a command-line interface (CLI) for quick conversions without writing a Python script.

```sh
python -m pymap3d geodetic2ecef 40.0 -105.0 1600.0
```

> -1266643.1360426995 -4727176.53876973 4079014.032375875

## PyMap3D functions

Functions include the
following, where the source coordinate system (before the "2") is
converted to the desired coordinate system:

```sh
python3 -m pymap3d -h
```

```
usage: pymap3d [-h]
               [--ellipsoid {maupertuis,plessis,everest1830,everest1830m,everest1967,airy,bessel,clarke1866,clarke1878,clarke1860,helmert,hayford,international1924,krassovsky1940,wgs66,australian,international1967,grs67,sa1969,wgs72,grs80,wgs84,wgs84_mean,iers1989,pz90.11,iers2003,gsk2011,mercury,venus,moon,mars,jupiter,io,saturn,uranus,neptune,pluto,jupyter}]
               [--output {text,json}] [--precision PRECISION]
               {aer2dca,aer2ecef,aer2eci,aer2enu,aer2geodetic,aer2ned,authalic2geodetic,azel2radec,conformal2geodetic,datetime2sidereal,dca2aer,dca2ecef,dca2enu,dca2geodetic,dca2ned,ecef2aer,ecef2dca,ecef2eci,ecef2enu,ecef2enuv,ecef2geodetic,ecef2ned,ecef2nedv,ecef2nvector,eci2aer,eci2ecef,eci2geodetic,enu2aer,enu2dca,enu2ecef,enu2ecefv,enu2geodetic,enu2uvw,geoc2geod,geocentric2geodetic,geocentric_radius,geod2geoc,geodetic2aer,geodetic2authalic,geodetic2conformal,geodetic2dca,geodetic2ecef,geodetic2eci,geodetic2enu,geodetic2geocentric,geodetic2isometric,geodetic2ned,geodetic2nvector,geodetic2parametric,geodetic2rectifying,geodetic2spherical,greenwichsrt,isometric2geodetic,meridian,ned2aer,ned2dca,ned2ecef,ned2geodetic,nvector2ecef,nvector2geodetic,parallel,parametric2geodetic,radec2azel,rectifying2geodetic,spherical2geodetic,str2dt,transverse,uvw2enu} ...

pymap3d CLI - Geographic coordinate conversions

positional arguments:
  {aer2dca,aer2ecef,aer2eci,aer2enu,aer2geodetic,aer2ned,authalic2geodetic,azel2radec,conformal2geodetic,datetime2sidereal,dca2aer,dca2ecef,dca2enu,dca2geodetic,dca2ned,ecef2aer,ecef2dca,ecef2eci,ecef2enu,ecef2enuv,ecef2geodetic,ecef2ned,ecef2nedv,ecef2nvector,eci2aer,eci2ecef,eci2geodetic,enu2aer,enu2dca,enu2ecef,enu2ecefv,enu2geodetic,enu2uvw,geoc2geod,geocentric2geodetic,geocentric_radius,geod2geoc,geodetic2aer,geodetic2authalic,geodetic2conformal,geodetic2dca,geodetic2ecef,geodetic2eci,geodetic2enu,geodetic2geocentric,geodetic2isometric,geodetic2ned,geodetic2nvector,geodetic2parametric,geodetic2rectifying,geodetic2spherical,greenwichsrt,isometric2geodetic,meridian,ned2aer,ned2dca,ned2ecef,ned2geodetic,nvector2ecef,nvector2geodetic,parallel,parametric2geodetic,radec2azel,rectifying2geodetic,spherical2geodetic,str2dt,transverse,uvw2enu}
                        Conversion command
    aer2dca             Converts AER (Azimuth, Elevation, Range) coordinates to DCA (Downrange,
                        Crossrange, Above).
    aer2ecef            converts target azimuth, elevation, range from observer at lat0,lon0,alt0
                        to ECEF coordinates.
    aer2eci             gives ECI of a point from an observer at az, el, slant range
    aer2enu             Azimuth, Elevation, Slant range to target to East, North, Up
    aer2geodetic        gives geodetic coordinates of a point with az, el, range
    aer2ned             converts azimuth, elevation, range to target from observer to North, East,
                        Down
    authalic2geodetic   converts from authalic latitude to geodetic latitude
    azel2radec          viewing angle (az, el) to sky coordinates (ra, dec)
    conformal2geodetic  converts from conformal latitude to geodetic latitude
    datetime2sidereal   Convert ``datetime`` to local sidereal time
    dca2aer             Converts DCA (Downrange, Crossrange, Above) coordinates to AER (Azimuth,
                        Elevation, Range).
    dca2ecef            Converts DCA (Downrange, Crossrange, Above) coordinates to ECEF (Earth-
                        Centered, Earth-Fixed) coordinates.
    dca2enu             Converts DCA (Downrange, Crossrange, Above) coordinates to ENU (East,
                        North, Up).
    dca2geodetic        Converts DCA (Downrange, Crossrange, Above) coordinates to geodetic
                        coordinates (latitude, longitude, altitude).
    dca2ned             Converts DCA (Downrange, Crossrange, Above) coordinates to NED (North,
                        East, Down).
    ecef2aer            compute azimuth, elevation and slant range from an Observer to a Point
                        with ECEF coordinates.
    ecef2dca            Converts ECEF (Earth-Centered, Earth-Fixed) coordinates to DCA (Downrange,
                        Crossrange, Above).
    ecef2eci            Point => Point ECEF => ECI
    ecef2enu            from observer to target, ECEF => ENU
    ecef2enuv           VECTOR from observer to target ECEF => ENU
    ecef2geodetic       convert ECEF (meters) to geodetic coordinates
    ecef2ned            Convert ECEF x,y,z to North, East, Down
    ecef2nedv           for VECTOR between two points
    ecef2nvector        Convert ECEF coordinates to an n-vector.
    eci2aer             takes Earth Centered Inertial x,y,z ECI coordinates of point and gives az,
                        el, slant range from Observer
    eci2ecef            Observer => Point ECI => ECEF
    eci2geodetic        convert Earth Centered Internal ECI to geodetic coordinates
    enu2aer             ENU to Azimuth, Elevation, Range
    enu2dca             Converts ENU (East, North, Up) coordinates to DCA (Downrange, Crossrange,
                        Above).
    enu2ecef            ENU to ECEF
    enu2ecefv           VECTOR from observer to target ENU => ECEF
    enu2geodetic        East, North, Up to target to geodetic coordinates
    enu2uvw             Parameters
    geoc2geod           convert geocentric latitude to geodetic latitude, consider mean sea level
                        altitude
    geocentric2geodetic
                        converts from geocentric latitude to geodetic latitude
    geocentric_radius   compute geocentric radius at geodetic latitude
    geod2geoc           convert geodetic latitude to geocentric latitude on spheroid surface
    geodetic2aer        gives azimuth, elevation and slant range from an Observer to a Point with
                        geodetic coordinates.
    geodetic2authalic   converts from geodetic latitude to authalic latitude
    geodetic2conformal  converts from geodetic latitude to conformal latitude
    geodetic2dca        Converts geodetic coordinates (latitude, longitude, altitude) to DCA
                        (Downrange, Crossrange, Above) coordinates.
    geodetic2ecef       point transformation from Geodetic of specified ellipsoid (default WGS-84)
                        to ECEF
    geodetic2eci        convert geodetic coordinates to Earth Centered Internal ECI
    geodetic2enu        Parameters
    geodetic2geocentric
                        convert geodetic latitude to geocentric latitude on spheroid surface
    geodetic2isometric  computes isometric latitude on an ellipsoid
    geodetic2ned        convert latitude, longitude, altitude of target to North, East, Down from
                        observer
    geodetic2nvector    Convert geodetic coordinates (latitude, longitude) to an n-vector.
    geodetic2parametric
                        converts from geodetic latitude to parametric latitude
    geodetic2rectifying
                        converts from geodetic latitude to rectifying latitude
    geodetic2spherical  point transformation from Geodetic of specified ellipsoid (default WGS-84)
    greenwichsrt        Convert Julian time to sidereal time
    isometric2geodetic  converts from isometric latitude to geodetic latitude
    meridian            computes the meridional radius of curvature for the ellipsoid
    ned2aer             converts North, East, Down to azimuth, elevation, range
    ned2dca             Converts NED (North, East, Down) coordinates to DCA (Downrange,
                        Crossrange, Above).
    ned2ecef            North, East, Down to target ECEF coordinates
    ned2geodetic        Converts North, East, Down to target latitude, longitude, altitude
    nvector2ecef        Convert an n-vector to ECEF coordinates.
    nvector2geodetic    Convert an n-vector back to geodetic coordinates (latitude, longitude).
    parallel            computes the radius of the small circle encompassing the globe at the
                        specified latitude
    parametric2geodetic
                        converts from parametric latitude to geodetic latitude
    radec2azel          sky coordinates (ra, dec) to viewing angle (az, el)
    rectifying2geodetic
                        converts from rectifying latitude to geodetic latitude
    spherical2geodetic  point transformation from geocentric spherical of specified ellipsoid
    str2dt              Converts times in string or list of strings to datetime(s)
    transverse          computes the radius of the curve formed by a plane
    uvw2enu             Parameters

options:
  -h, --help            show this help message and exit
  --ellipsoid {maupertuis,plessis,everest1830,everest1830m,everest1967,airy,bessel,clarke1866,clarke1878,clarke1860,helmert,hayford,international1924,krassovsky1940,wgs66,australian,international1967,grs67,sa1969,wgs72,grs80,wgs84,wgs84_mean,iers1989,pz90.11,iers2003,gsk2011,mercury,venus,moon,mars,jupiter,io,saturn,uranus,neptune,pluto,jupyter}
                        Ellipsoid model to use (default: wgs84)
  --output {text,json}  Output format (default: text)
  --precision PRECISION
                        Decimal places for floating-point text output (default: None)
```

Vincenty functions "vincenty.vreckon" and "vincenty.vdist" are accessed like:

```python
import pymap3d.vincenty as pmv

lat2, lon2 = pmv.vreckon(lat1, lon1, ground_range_m, azimuth_deg)
dist_m, azimuth_deg = pmv.vdist(lat1, lon1, lat2, lon2)
```

Additional functions:

* loxodrome_inverse: rhumb line distance and azimuth between ellipsoid points (lat,lon)  akin to Matlab `distance('rh', ...)` and `azimuth('rh', ...)`
* loxodrome_direct
* geodetic latitude transforms to/from: parametric, authalic, isometric, and more in pymap3d.latitude

Abbreviations:

* [AER: Azimuth, Elevation, Range](https://en.wikipedia.org/wiki/Spherical_coordinate_system)
* [ECEF: Earth-centered, Earth-fixed](https://en.wikipedia.org/wiki/ECEF)
* [ECI: Earth-centered Inertial using IERS](https://www.iers.org/IERS/EN/Home/home_node.html) via `astropy`
* [ENU: East North Up](https://en.wikipedia.org/wiki/Axes_conventions#Ground_reference_frames:_ENU_and_NED)
* [NED: North East Down](https://en.wikipedia.org/wiki/North_east_down)
* [radec: right ascension, declination](https://en.wikipedia.org/wiki/Right_ascension)

### Ellipsoid

Numerous functions in pymap3d use an ellipsoid model.
The default is WGS84 Ellipsoid.
Numerous other ellipsoids are available in pymap3d.Ellipsoid.

Print available ellipsoid models:

```python
import pymap3d as pm

print(pm.Ellipsoid.models)
```

Specify GRS80 ellipsoid:

```python
import pymap3d as pm

ell = pm.Ellipsoid.from_name('grs80')
```

### array vs scalar

Use of pymap3d on embedded systems or other streaming data applications often deal with scalar position data.
These data are handled efficiently with the Python math stdlib module.
Vector data can be handled via list comprehension.

Those needing multidimensional data with SIMD and other Numpy and/or PyPy accelerated performance can do so automatically by installing Numpy.
pymap3d seamlessly falls back to Python's math module if Numpy isn't present.
To keep the code clean, only scalar data can be used without Numpy.
As noted above, use list comprehension if you need vector data without Numpy.

### Caveats

* Atmospheric effects neglected in all functions not invoking AstroPy.
  Would need to update code to add these input parameters (just start a GitHub Issue to request).
* Planetary perturbations and nutation etc. not fully considered.

## Compare to Matlab Mapping and Aerospace Toolbox

The tests in files tests/test_matlab*.py selected by

```sh
pytest -k matlab
# run from pymap3d/ top-level directory
```

use
[Matlab Engine for Python](https://www.mathworks.com/help/matlab/matlab_external/install-the-matlab-engine-for-python.html)
to compare Python PyMap3D output with Matlab output using Matlab functions.

```sh
python -m pip install matlabengine
```

## Notes

As compared to [PyProj](https://github.com/jswhit/pyproj):

* PyMap3D does not require anything beyond pure Python for most transforms
* Astronomical conversions are done using (optional) AstroPy for established accuracy
* PyMap3D API is similar to Matlab Mapping Toolbox, while PyProj's interface is quite distinct
* PyMap3D intrinsically handles local coordinate systems such as ENU,
  while PyProj ENU requires some [additional effort](https://github.com/jswhit/pyproj/issues/105).
* PyProj is oriented towards points on the planet surface, while PyMap3D handles points on or above the planet surface equally well, particularly important for airborne vehicles and remote sensing.

### AstroPy.Units.Quantity

At this time,
[AstroPy.Units.Quantity](http://docs.astropy.org/en/stable/units/)
is not supported.
Let us know if this is of interest.
Impacts on performance would have to be considered before making Quantity a first-class citizen.
For now, you can workaround by passing in the `.value` of the variable.
