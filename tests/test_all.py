#!/usr/bin/env python
# Copyright (c) 2014-2018 Michael Hirsch, Ph.D.
#
# Redistribution and use in source and binary forms, with or without modification, are permitted provided that the following conditions are met:
# 1. Redistributions of source code must retain the above copyright notice, this list of conditions and the following disclaimer.
# 2. Redistributions in binary form must reproduce the above copyright notice, this list of conditions and the following disclaimer in the documentation and/or other materials provided with the distribution.
# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

"""
runs tests
"""
import logging
import unittest
from datetime import datetime
try:
    import numpy, astropy
    from numpy import radians
    from pymap3d.vincenty import vreckon,vdist
except ImportError:
    numpy = astropy = None
    from math import radians

try: # for validation
    import pyproj
except ImportError:
    pyproj = None

from pytz import UTC
#
import pymap3d as pm
from pymap3d.haversine import anglesep,anglesep_meeus
from pymap3d.datetime2hourangle import datetime2sidereal
from pymap3d.vallado import vazel2radec, vradec2azel
from pymap3d.timeconv import str2dt

t0 = '2014-04-06T08:00:00Z'
lat, lon = (65, -148)
ra, dec = (166.5032081149338,55.000011165405752)
ha = 45.482789587392013
azi, eli = (180.1,80)

lla0 = 42, -82, 200
rlla0 = (radians(lla0[0]), radians(lla0[1]), lla0[2])

aer0 = (33, 70, 1000)
raer0 = (radians(aer0[0]), radians(aer0[1]), aer0[2])

# %% outcomes from matlab
xyz0 = (660.6753e3, -4700.9487e3, 4245.738e3) # geodetic2ecef

lla1 = (42.002582, -81.997752, 1.1397018e3) #aer2geodetic
rlla1 = (radians(lla1[0]), radians(lla1[1]), lla1[2])

axyz0 = 660930.2, -4701424, 4246579.6 #aer2ecef

enu0 = (186.277521, 286.842228, 939.692621) #aer2enu
ned0 = (enu0[1],enu0[0],-enu0[2])

# vector
vx,vy,vz = (5,3,2)
ve,vn,vu =(5.368859646588048, 3.008520763668120, -0.352347711524077)



class Pure(unittest.TestCase):

    def test_str2dt(self):

        assert str2dt(datetime(2014,4,6,8, tzinfo=UTC)) == datetime(2014,4,6,8, tzinfo=UTC) # passthrough
        assert str2dt('2014-04-06T08:00:00Z') == datetime(2014,4,6,8, tzinfo=UTC)
        ti = [str2dt('2014-04-06T08:00:00Z'), str2dt('2014-04-06T08:01:02Z')]
        to = [datetime(2014,4,6,8, tzinfo=UTC),  datetime(2014,4,6,8,1,2, tzinfo=UTC)]
        assert  ti == to   # even though ti is numpy array of datetime and to is list of datetime

# %%

    def test_datetime2sidereal(self):
        sra = 2.90658
        # http://www.jgiesen.de/astro/astroJS/siderealClock/
        try:
            assert_allclose(datetime2sidereal(t0, radians(lon), False), sra, rtol=1e-5)

            assert_allclose(datetime2sidereal([t0], radians(lon), False), [sra], rtol=1e-5)
        except ImportError:
            pass

        assert_allclose(datetime2sidereal(t0, radians(lon), True), sra, rtol=1e-5)

        assert_allclose(datetime2sidereal([t0], radians(lon), True), [sra], rtol=1e-5)



    def test_haversine(self):
        try:
            assert_allclose(anglesep(35,23, 84,20), ha)
        except ImportError:
            pass
        #%% compare with astropy
        assert_allclose(anglesep_meeus(35,23, 84,20), ha)


    def test_geodetic(self):
        if pyproj:
            ecef = pyproj.Proj(proj='geocent', ellps='WGS84', datum='WGS84')
            lla = pyproj.Proj(proj='latlong', ellps='WGS84', datum='WGS84')

        xyz1 = pm.geodetic2ecef(*lla0)

        assert_allclose(pm.geodetic2ecef(*rlla0, deg=False), xyz1, err_msg='geodetic2ecef: rad')
        assert_allclose(xyz1, xyz0, err_msg='geodetic2ecef: deg')

        assert_allclose(pm.ecef2geodetic(*xyz1), lla0, err_msg='ecef2geodetic: deg')
        assert_allclose(pm.ecef2geodetic(*xyz1, deg=False), rlla0, err_msg='ecef2geodetic: rad')

        if pyproj:
            assert_allclose(pyproj.transform(lla,ecef,lla0[1],lla0[0],lla0[2]), xyz1)
            assert_allclose(pyproj.transform(ecef, lla, *xyz1), 
                            (lla0[1],lla0[0],lla0[2]))


        lla2 = pm.aer2geodetic(*aer0, *lla0)
        rlla2 = pm.aer2geodetic(*raer0, *rlla0, deg=False)

        assert_allclose(lla2, lla1, err_msg='aer2geodetic: deg')
        assert_allclose(rlla2, rlla1, err_msg='aer2geodetic:rad')

        assert_allclose(pm.geodetic2aer(*lla2, *lla0), aer0, err_msg= 'geodetic2aer: deg')
        assert_allclose(pm.geodetic2aer(*rlla2, *rlla0, deg=False), raer0, err_msg= 'geodetic2aer: rad')
        
# %% aer-ecef
        xyz2 = pm.aer2ecef(*aer0,*lla0)

        assert_allclose(pm.aer2ecef(*raer0, *rlla0,deg=False),
                                    axyz0, err_msg='aer2ecef:rad')

        assert_allclose(xyz2, axyz0, err_msg='aer2ecef: deg')

        assert_allclose(pm.ecef2aer(*xyz2, *lla0), aer0, err_msg='ecef2aer:deg')
        assert_allclose(pm.ecef2aer(*xyz2, *rlla0, deg=False), raer0, err_msg='ecef2aer:rad')
# %% aer-enu
        enu1 = pm.aer2enu(*aer0)
        ned1 = (enu1[1],enu1[0],-enu1[2])
        
        assert_allclose(enu1,enu0, err_msg='aer2enu: deg')
        assert_allclose(pm.aer2enu(*raer0, deg=False), enu0, err_msg='aer2enu: rad')

        assert_allclose(pm.aer2ned(*aer0), ned0, err_msg='aer2ned')

        assert_allclose(pm.enu2aer(*enu1), aer0, err_msg='enu2aer: deg')
        assert_allclose(pm.enu2aer(*enu1, deg=False), raer0, err_msg='enu2aer: rad')

        assert_allclose(pm.ned2aer(*ned1), aer0, err_msg='ned2aer')

# %% enu-ecef
        assert_allclose(pm.enu2ecef(*enu1, *lla0), xyz2, err_msg='enu2ecef: deg')
        assert_allclose(pm.enu2ecef(*enu1, *rlla0, deg=False), xyz2, err_msg='enu2ecef: rad')

        assert_allclose(pm.ecef2enu(*xyz2, *lla0), enu1, err_msg='ecef2enu:deg')
        assert_allclose(pm.ecef2enu(*xyz2, *rlla0, deg=False), enu1, err_msg='ecef2enu:rad')

        assert_allclose(pm.ecef2ned(*xyz2, *lla0), ned1,
                        err_msg='ecef2ned')

        assert_allclose(pm.ned2ecef(*ned1, *lla0), xyz2,
                        err_msg='ned2ecef')
# %%
        assert_allclose(pm.ecef2enuv(vx,vy,vz, *lla0[:2]), (ve,vn,vu))


        assert_allclose(pm.ecef2nedv(vx,vy,vz, *lla0[:2]), (vn,ve,-vu))

# %%
        enu3 = pm.geodetic2enu(*lla2, *lla0)
        ned3 = (enu3[1],enu3[0],-enu3[2])

        assert_allclose(pm.geodetic2ned(*lla2, *lla0), ned3,
                        err_msg='geodetic2ned: deg')

        assert_allclose(pm.enu2geodetic(*enu3, *lla0), lla2,
                        err_msg='enu2geodetic')

        assert_allclose(pm.ned2geodetic(*ned3, *lla0), lla2,
                        err_msg='ned2geodetic')

class Numpy(unittest.TestCase):

    def test_vincenty(self):
        if numpy is None:
            logging.warning('Vincenty not tested')
            return

        az = 38; sr = 3e3
        lat2,lon2,a21 = vreckon(10,20,sr,az)
        assert_allclose((lat2,lon2,a21),
                            (10.02137267,20.016847,218.0029286))
        if pyproj:
            p4lon,p4lat,p4a21 = pyproj.Geod(ellps='WGS84').fwd(lon2,lat2,az,sr)
            assert_allclose((p4lon,p4lat,p4a21 % 360.), (lon2,lat2, a21), rtol=0.0025)


        assert_allclose(vdist(10,20,lat2,lon2), (sr,az,a21))
        if pyproj:
            p4az,p4a21,p4sr = pyproj.Geod(ellps='WGS84').inv(20,10,lon2,lat2)
            assert_allclose((p4az,p4a21 % 360.,p4sr), (az,a21,sr))


    def test_azel2radec(self):
        R,D = pm.azel2radec(azi, eli, lat, lon, t0)
        assert_allclose(R, ra, rtol=1e-2)
        assert_allclose(D, dec, rtol=1e-2)

        Rv, Dv = vazel2radec(azi, eli, lat, lon, t0)
        assert_allclose(Rv, ra)
        assert_allclose(Dv, dec)


    def test_radec2azel(self):
        if numpy is None:
            logging.warning('RA DEC not tested')
            return
        azapy, elapy = pm.radec2azel(ra,dec,  lat, lon, t0)
        assert_allclose(azapy, azi, rtol=1e-2)
        assert_allclose(elapy, eli, rtol=1e-2)

        azvallado, elvallado = vradec2azel(ra, dec, lat, lon, t0)
        assert_allclose(azvallado, azi, rtol=1e-2)
        assert_allclose(elvallado, eli, rtol=1e-2)



    def test_eci(self):
        if numpy is None or astropy is None:
            logging.warning('ECI not tested')
            return
            
        teci = (-3.977913815668146e6,-2.582332196263046e6,4.250818828152067e6)
        t = datetime(2013,1,15,12,0,5,tzinfo=UTC)
        lla = numpy.asarray(pm.eci2geodetic(teci,t)).squeeze()
        assert_allclose(lla, lla0,rtol=0.2)

        assert_allclose(pm.eci2ecef(teci,t).squeeze(),
                    [649012.04640917,-4697980.55129606,4250818.82815207])

        assert_allclose(pm.ecef2eci([649012.04640917,-4697980.55129606,4250818.82815207],t).squeeze(),
                        teci)

        assert_allclose(numpy.asarray(pm.eci2aer(teci,42,-100,0,t)).squeeze(),
                        [83.73050,-6.614478,1.473510e6])

#%%

def assert_allclose(actual, desired, rtol=1e-7, atol=0,err_msg=''):
    if numpy:
        numpy.testing.assert_allclose(actual, desired, rtol, atol)
    else:
        assert pm.allclose(actual, desired, rtol, atol)


if __name__ == '__main__':
    unittest.main()
