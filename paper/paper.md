---
title: 'PyMap3D: 3-D coordinate conversions for terrestrial and geospace environments'
tags:
authors:
  - name: Michael Hirsch
    orcid: 0000-0002-1637-6526
    affiliation: "1, 2"
affiliations:
 - name: Boston University ECE Dept.
   index: 1
 - name: SciVision, Inc.
   index: 2
date: 29 January 2018
bibliography: paper.bib
---

# Summary

PyMap3D [@pymap3d] is a Python coordinate transformation program that converts between geospace coordinate systems.
A subset of the Python functionality using syntax compatible with the $1000 Matlab Mapping Toolbox is provided for Matlab and GNU Octave users in the ``matlab/`` directory.
Increased performance and accuracy is optionally available for certain functions for those using AstroPy, but it is not required.
Since Numpy is used throughout PyMap3D, many of the functions operate with speed comparable to compiled languages.
Other functions that are iterative could possibly be sped up with modules such as Cython or Numba.

PyMap3D is targeted for users needing conversions between coordinate systems for observation platforms near Earth's surface, whether underwater, ground-based or space-based platforms.
This includes rocket launches, orbiting spacecrafts, UAVs, cameras, radars and many more. 
By adding ellipsoid parameters, it could be readily be used for other planets as well.
The coordinate systems included are:
* ECEF (Earth centered, Earth fixed)
* ENU (East, North, Up)
* NED (North, East, Down)
* ECI (Earth Centered Inertial)
* Geodetic (Latitude, Longitude, Altitude)
* Horizontal Celestial (Alt-Az or Az-El)
* Equatorial Celestial (Right Ascension, Declination)

Additionally, Vincenty [@vincenty, @veness] geodesic distances and direction are computed.

PyMap3D has already seen usage in projects including
* [EU ECSEL project 662107 SWARMs](http://swarms.eu/)
* Rafael Defense Systems DataHack 2017
* HERA radiotelescope
* Mahali (NSF Grant: AGS-1343967)
* Solar Eclipse network (NSF Grant: AGS-1743832)
* High Speed Auroral Tomography (NSF Grant: AGS-1237376) [@7368896]

## Other Programs

Other Python geodesy programs include:

* [PyGeodesy](https://github.com/mrJean1/PyGeodesy) MIT license
* [PyProj](https://github.com/jswhit/pyproj) ISC license

These programs are targeted for geodesy experts, and require additional packages beyond Python that may not be readily accessible to users.
Further, these programs do not include all the functions of PyMap3D, and do not have the straightforward function-based API of PyMap3D.




# References
