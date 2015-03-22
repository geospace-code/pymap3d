
[![Build Status](https://travis-ci.org/scienceopen/python-mapping.svg)](https://travis-ci.org/scienceopen/python-mapping)
python-mapping
==============

Python coordinate conversions simply ported from Matlab/Octave functions
Credit goes to original .m file authors.

Prereq:
```
conda install --file requirements.txt
```
or
```
pip install -r requirements.txt
```

Installation:
-------------
```
pip install https://github.com/scienceopen/python-mapping/archive/master.zip
```

Consider using python-geopy as well.

Popular mapping toolbox functions ported to Python include:
```
aer2ecef
aer2enu
aer2geodetic
aer2ned
ecef2aer
ecef2enu
ecef2enuv
ecef2geodetic
ecef2ned
ecef2nedv
ecef2eci
eci2ecef
enu2aer
enu2ecef
enu2ecefv
enu2geodetic
geodetic2aer
geodetic2ecef
geodetic2enu
geodetic2ned
ned2aer
ned2ecef
ned2ecefv
ned2geodetic 
vreckon
```
for converting right ascension and declination to azimuth and elevation, please see the function radec2azel inside

https://github.com/scienceopen/astrometry/

or simply use the functionality of AstroPy 1.0+ to do radec->azel conversion 

http://astropy.readthedocs.org/en/v1.0/whatsnew/1.0.html#support-for-alt-az-coordinates
