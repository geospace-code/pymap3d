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
from pytz import UTC
#
import pymap3d as pm
from pymap3d.haversine import anglesep,anglesep_meeus
from pymap3d.datetime2hourangle import datetime2sidereal
from pymap3d.vallado import vazel2radec, vradec2azel
from pymap3d.timeconv import str2dt

t = '2014-04-06T08:00:00Z'
lat, lon = (65, -148)
ra, dec = (166.5032081149338,55.000011165405752)
ha = 45.482789587392013
azi, eli = (180.1,80)

tlat,tlon,talt = 42, -82, 200
taz,tel,tsrange = 33, 70, 1000
# %% outcomes from matlab
x0, y0, z0 = 660.6753e3, -4700.9487e3, 4245.738e3 # geodetic2ecef
lat1, lon1, alt1 = 42.002582, -81.997752, 1.1397018e3 #aer2geodetic
a2x, a2y, a2z = 660930.2, -4701424, 4246579.6 #aer2ecef
e0, n0, u0 = 186.277521, 286.842228, 939.692621 #aer2enu

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
            sdrapy = datetime2sidereal(t, radians(lon), False)
            assert_allclose(sdrapy, sra, rtol=1e-5)

            sdrapy = datetime2sidereal([t], radians(lon), False)
            assert_allclose(sdrapy, [sra], rtol=1e-5)
        except ImportError:
            pass

        sdrvallado = datetime2sidereal(t, radians(lon), True)
        assert_allclose(sdrvallado, sra, rtol=1e-5)

        sdrvallado = datetime2sidereal([t], radians(lon), True)
        assert_allclose(sdrvallado, [sra], rtol=1e-5)



    def test_haversine(self):
        try:
            assert_allclose(anglesep(35,23, 84,20), ha)
        except ImportError:
            pass
        #%% compare with astropy
        assert_allclose(anglesep_meeus(35,23, 84,20), ha)


    def test_geodetic(self):
        x1,y1,z1 = pm.geodetic2ecef(tlat,tlon,talt)

        assert_allclose(pm.geodetic2ecef(radians(tlat),radians(tlon),talt,deg=False),
                        (x1,y1,z1))

        assert_allclose((x1,y1,z1), (x0,y0,z0),
                        err_msg='geodetic2ecef:')

        assert_allclose(pm.ecef2geodetic(x1,y1,z1), (tlat,tlon,talt),
                    err_msg='ecef2geodetic:')


        lat2,lon2,alt2 = pm.aer2geodetic(taz,tel,tsrange,tlat,tlon,talt)

        assert_allclose((lat2,lon2,alt2), (lat1,lon1,alt1),
                   err_msg='aer2geodetic')

        assert_allclose(pm.geodetic2aer(lat2,lon2,alt2,tlat,tlon,talt),
                        (taz,tel,tsrange),
                         err_msg= 'geodetic2aer')


        x2,y2,z2 = pm.aer2ecef(taz,tel,tsrange,tlat,tlon,talt)

        assert_allclose(pm.aer2ecef(radians(taz),radians(tel),tsrange,radians(tlat),radians(tlon),talt,deg=False),
                        (a2x,a2y,a2z))

        assert_allclose((x2,y2,z2), (a2x,a2y,a2z),
                         err_msg='aer2ecef')

        assert_allclose(pm.ecef2aer(x2, y2, z2, tlat, tlon,talt), (taz,tel,tsrange),
                        err_msg='ecef2aer')


        e1,n1,u1 = pm.aer2enu(taz,tel,tsrange)

        assert_allclose((e1,n1,u1),(e0,n0,u0),
                        err_msg='aer2enu')

        assert_allclose(pm.aer2ned(taz,tel,tsrange),(n0,e0,-u0),
                        err_msg='aer2ned')

        assert_allclose(pm.enu2aer(e1,n1,u1), (taz,tel,tsrange),
                        err_msg='enu2aer')

        assert_allclose(pm.ned2aer(n1,e1,-u1), (taz,tel,tsrange),
                        err_msg='ned2aer')


        assert_allclose(pm.enu2ecef(e1,n1,u1,tlat,tlon,talt),(x2, y2, z2),
                        err_msg='enu2ecef')

        assert_allclose(pm.ecef2enu(x2,y2,z2, tlat, tlon, talt),(e1,n1,u1),
                        err_msg='ecef2enu')

        assert_allclose(pm.ecef2ned(x2,y2,z2, tlat, tlon, talt),(n1,e1,-u1),
                        err_msg='ecef2ned')

        assert_allclose(pm.ned2ecef(n1,e1,-u1,tlat,tlon,talt),(x2,y2,z2),
                        err_msg='ned2ecef')
    # %%
        assert_allclose(pm.ecef2enuv(vx,vy,vz,tlat,tlon), (ve,vn,vu))


        assert_allclose(pm.ecef2nedv(vx,vy,vz,tlat,tlon), (vn,ve,-vu))

    #%%
        e3,n3,u3 = pm.geodetic2enu(lat2,lon2,alt2,tlat,tlon,talt)

        assert_allclose(pm.geodetic2ned(lat2,lon2,alt2,tlat,tlon,talt),(n3,e3,-u3))

        assert_allclose(pm.enu2geodetic(e3,n3,u3,tlat,tlon,talt),(lat2,lon2,alt2),
                        err_msg='enu2geodetic')

        assert_allclose(pm.ned2geodetic(n3,e3,-u3,tlat,tlon,talt),(lat2,lon2,alt2),
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

        assert_allclose(vdist(10,20,lat2,lon2),(sr,az,a21))


    def test_azel2radec(self):
        R,D = pm.azel2radec(azi, eli, lat, lon, t)
        assert_allclose(R, ra, rtol=1e-2)
        assert_allclose(D, dec, rtol=1e-2)

        Rv, Dv = vazel2radec(azi, eli, lat, lon, t)
        assert_allclose(Rv, ra)
        assert_allclose(Dv, dec)


    def test_radec2azel(self):
        if numpy is None:
            logging.warning('RA DEC not tested')
            return
        azapy, elapy = pm.radec2azel(ra,dec,  lat, lon, t)
        assert_allclose(azapy, azi, rtol=1e-2)
        assert_allclose(elapy, eli, rtol=1e-2)

        azvallado, elvallado = vradec2azel(ra, dec, lat, lon, t)
        assert_allclose(azvallado, azi, rtol=1e-2)
        assert_allclose(elvallado, eli, rtol=1e-2)



    def test_eci(self):
        if numpy is None or astropy is None:
            logging.warning('ECI not tested')
            return
        tlla = (tlat, tlon, talt)
        teci = (-3.977913815668146e6,-2.582332196263046e6,4.250818828152067e6)
        t = datetime(2013,1,15,12,0,5,tzinfo=UTC)
        lla = numpy.asarray(pm.eci2geodetic(teci,t)).squeeze()
        assert_allclose(lla,tlla,rtol=0.2)

        assert_allclose(pm.eci2ecef(teci,t).squeeze(),
                    [649012.04640917,-4697980.55129606,4250818.82815207])

        assert_allclose(pm.ecef2eci([649012.04640917,-4697980.55129606,4250818.82815207],t).squeeze(),
                        teci)

        assert_allclose(numpy.asarray(pm.eci2aer(teci,42,-100,0,t)).squeeze(),
                        [83.73050,-6.614478,1.473510e6])

#%%

def assert_allclose(actual, desired, rtol=1e-7, atol=0,err_msg=''):
    assert pm.allclose(actual, desired, rtol, atol)


if __name__ == '__main__':
    unittest.main()
