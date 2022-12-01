# Python 3-D coordinate conversions

[![image](https://zenodo.org/badge/DOI/10.5281/zenodo.213676.svg)](https://doi.org/10.5281/zenodo.213676)
[![image](http://joss.theoj.org/papers/10.21105/joss.00580/status.svg)](https://doi.org/10.21105/joss.00580)
[![codecov](https://codecov.io/gh/geospace-code/pymap3d/branch/main/graph/badge.svg?token=DFWBW6TKNr)](https://codecov.io/gh/geospace-code/pymap3d)
![Actions Status](https://github.com/geospace-code/pymap3d/workflows/ci/badge.svg)
![Actions Status](https://github.com/geospace-code/pymap3d/workflows/ci_stdlib_only/badge.svg)
[![image](https://img.shields.io/pypi/pyversions/pymap3d.svg)](https://pypi.python.org/pypi/pymap3d)
[![PyPi Download stats](http://pepy.tech/badge/pymap3d)](http://pepy.tech/project/pymap3d)

Pure Python (no prerequistes beyond Python itself) 3-D geographic coordinate conversions and geodesy.
API similar to popular $1000 Matlab Mapping Toolbox routines for Python
PyMap3D is intended for non-interactive use on massively parallel (HPC) and embedded systems.

[API docs](https://geospace-code.github.io/pymap3d/)

Thanks to our [contributors](./.github/contributors.md).

## Similar toolboxes in other code languages

* [Matlab, GNU Octave](https://github.com/geospace-code/matmap3d)
* [Fortran](https://github.com/geospace-code/maptran3d)
* [Rust](https://github.com/gberrante/map_3d)

## Prerequisites

Pymap3d is compatible with Python &ge; 3.7 including PyPy.
Numpy and AstroPy are optional; algorithms from Vallado and Meeus are used if AstroPy is not present.

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
obslla = (obs_lat,obs_lon,obs_alt)

lla = pm.aer2geodetic(*aer,*obslla)
```

where tuple `lla` is comprised of scalar or N-D arrays `(lat,lon,alt)`.

Example scripts are in the [examples](./Examples) directory.

Native Python float is typically [64 bit](https://docs.python.org/3/library/stdtypes.html#typesnumeric).
Numpy can select real precision bits: 32, 64, 128, etc.

### Functions

Popular mapping toolbox functions ported to Python include the
following, where the source coordinate system (before the "2") is
converted to the desired coordinate system:

```
aer2ecef  aer2enu  aer2geodetic  aer2ned
ecef2aer  ecef2enu  ecef2enuv  ecef2geodetic  ecef2ned  ecef2nedv
ecef2eci  eci2ecef eci2aer aer2eci geodetic2eci eci2geodetic
enu2aer  enu2ecef   enu2geodetic
geodetic2aer  geodetic2ecef  geodetic2enu  geodetic2ned
ned2aer  ned2ecef   ned2geodetic
azel2radec radec2azel
lookAtSpheroid
track2 departure meanm
rcurve rsphere
geod2geoc geoc2geod
geodetic2spherical spherical2geodetic
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
