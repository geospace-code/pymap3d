#!/usr/bin/env python
"""
runs tests
"""
import subprocess
from pathlib import Path
from datetime import datetime
from numpy import asarray,radians
from numpy.testing import assert_allclose, assert_almost_equal,run_module_suite
from pytz import UTC
#
import pymap3d as pm
from pymap3d.haversine import angledist,angledist_meeus
from pymap3d.vincenty import vreckon,vdist
from pymap3d.datetime2hourangle import datetime2sidereal
from pymap3d.vallado import vazel2radec, vradec2azel
from pymap3d.timeconv import str2dt
#
path = Path(__file__).parents[1]

def test_bsr():
    subprocess.check_call(['octave-cli','-q','Test.m'], cwd=path/'tests')


def test_str2dt():

    assert str2dt('2014-04-06T08:00:00Z') == datetime(2014,4,6,8, tzinfo=UTC)
    ti = [str2dt('2014-04-06T08:00:00Z'), str2dt('2014-04-06T08:01:02Z')]
    to = [datetime(2014,4,6,8, tzinfo=UTC),  datetime(2014,4,6,8,1,2, tzinfo=UTC)]
    assert  ti == to   # even though ti is numpy array of datetime and to is list of datetime

# %%
t = '2014-04-06T08:00:00Z'
rdlat, rdlon = (65, -148)
ra, dec = (166.5032081149338,55.000011165405752)
ha = 45.482789587392013
azi, eli = (180.1,80)

def test_datetime2sidereal():
    sdrapy = datetime2sidereal(t, radians(rdlon), False)
    assert_allclose(sdrapy, 2.9065780550600806, rtol=1e-5)

    sdrvallado = datetime2sidereal(t, radians(rdlon), True)
    assert_allclose(sdrvallado, 2.9065780550600806, rtol=1e-5)


def test_azel2radec():
    R,D = pm.azel2radec(azi, eli, rdlat, rdlon, t)
    assert_allclose(R, ra, rtol=1e-2)
    assert_allclose(D, dec, rtol=1e-2)

    Rv, Dv = vazel2radec(azi, eli, rdlat, rdlon, t)
    assert_allclose(Rv, ra)
    assert_allclose(Dv, dec)


def test_radec2azel():
    azapy, elapy = pm.radec2azel(ra,dec,  rdlat, rdlon, t)
    assert_allclose(azapy, azi, rtol=1e-2)
    assert_allclose(elapy, eli, rtol=1e-2)

    azvallado, elvallado = vradec2azel(ra, dec, rdlat, rdlon, t)
    assert_allclose(azvallado, azi, rtol=1e-2)
    assert_allclose(elvallado, eli, rtol=1e-2)


def test_haversine():
    assert_almost_equal(angledist(35,23, 84,20), ha)
    #%% compare with astropy
    assert_almost_equal(ha, angledist_meeus(35,23, 84,20))


def test_vreckon():
    lat2,lon2,a21 = vreckon(10,20,3000,38)
    assert_almost_equal(lat2,10.021372672660874)
    assert_almost_equal(lon2,20.016847098929979)
    assert_almost_equal(a21,218.0029285624942)
    #assert vreckon(91,0,0,0) is None

#%% coordconv3d
tlat,tlon,talt = 42, -82, 200
lat2, lon2, alt2 = 42.1, -81.9, 1300
taz,tel,tsrange = 33, 70, 1000
tx, ty, tz  =  (6.678411289903646e+05,
             -4.692496355102768e+06,
             4.254052899714093e+06)
te, tn, tu = (8.273771039503677e+03,
             1.111452002615149e+04,
             1.084939260985176e+03)
#outcomes from matlab
a2e, a2n, a2u = 186.277521, 286.842228, 939.692621 #aer2enu
a2x, a2y, a2z = 660930.192761, -4701424.222957, 4246579.604633 #aer2ecef
a2la, a2lo, a2a = 42.002582, -81.997752, 1139.701800 #aer2geodetic

g2az, g2el, g2rn = 36.664403, 4.477195, 13898.378892 #geodetic2aer
ec2az, ec2el, ec2rn =  36.664419, 0.351293, 13854.054270 #ecef2aer

g2x, g2y, g2z = 660675.251825, -4700948.683162, 4245737.662222
#geodeteic2ecef
e2e, e2n, e2u = 8272.476048, 11112.773942, 84.941624 #ecef2enu

ec2la, ec2lo, ec2a = 42.100000, -81.900000, 300.000000 #ecef2geodetic

#enu2ecef
e2x, e2y, e2z = (6.679456743004259e+05, -4.693230928738789e+06, 4.254723326333052e+06)

e2az, e2el, e2rn = (36.664402767128749, 4.477194667550686, 1.389837889201037e+04)

# vector
vx,vy,vz = (5,3,2)
ve,vn,vu =(5.368859646588048, 3.008520763668120, -0.352347711524077)


def test_geodetic():
    assert_allclose(pm.ecef2geodetic(tx,ty,tz),(ec2la,ec2lo,ec2a),
                rtol=0.01,
                err_msg='ecef2geodetic: ' + str(pm.ecef2geodetic(tx,ty,tz))
                )

    assert_allclose(pm.geodetic2aer(lat2,lon2,alt2,tlat,tlon,talt),
                    (g2az,g2el,g2rn),
                    rtol=0.05,
                    err_msg= 'geodetic2aer: {}'.format(
                    pm.geodetic2aer(lat2,lon2,alt2,tlat,tlon,talt)))

    assert_allclose(pm.geodetic2ecef(tlat,tlon,talt),(g2x,g2y,g2z),
                    rtol=0.01,
                    err_msg='geodetic2ecef: {}'.format(
                    pm.geodetic2ecef(tlat,tlon,talt)))

def test_ecefenu():
    assert_allclose(pm.aer2ecef(taz,tel,tsrange,tlat,tlon,talt), (a2x,a2y,a2z),
                    rtol=0.01,
                    err_msg='aer2ecef: {}'.format(
                    pm.aer2ecef(taz,tel,tsrange,tlat,tlon,talt)))

    assert_allclose(pm.aer2enu(taz,tel,tsrange),(a2e,a2n,a2u),
                    rtol=0.01,
                    err_msg='aer2enu: ' + str(pm.aer2enu(taz,tel,tsrange)))

    assert_allclose(pm.aer2ned(taz,tel,tsrange),(a2n,a2e,-a2u),
                    rtol=0.01,
                    err_msg='aer2ned: ' + str(pm.aer2ned(taz,tel,tsrange)))


    assert_allclose(pm.ecef2enu(tx,ty,tz, tlat, tlon, talt),(e2e,e2n,e2u),
                    rtol=0.01,
                    err_msg='ecef2enu: {}'.format(pm.ecef2enu(tx,ty,tz, tlat,
                    tlon, talt)))

    assert_allclose(pm.ecef2enuv(vx,vy,vz,tlat,tlon), (ve,vn,vu))

    assert_allclose(pm.ecef2ned(tx,ty,tz, tlat, tlon, talt),(e2n,e2e,-e2u),
                    rtol=0.01,
                    err_msg='ecef2ned: {}'.format(pm.ecef2enu(tx,ty,tz, tlat,
                    tlon, talt)))

    assert_allclose(pm.ecef2nedv(vx,vy,vz,tlat,tlon), (vn,ve,-vu))

    assert_allclose(pm.aer2geodetic(taz,tel,tsrange,tlat,tlon,talt),(a2la,a2lo,a2a),
                rtol=0.01,err_msg='aer2geodetic {}'.format(
                  pm.aer2geodetic(taz,tel,tsrange,tlat,tlon,talt)))

    assert_allclose(pm.ecef2aer(tx, ty, tz, tlat, tlon,talt),
                    (ec2az,ec2el,ec2rn),
                    rtol=0.01,
                    err_msg='ecef2aer {}'.format(pm.ecef2aer(a2x, a2y, a2z, tlat,
                    tlon, talt)))
#%%
    assert_allclose(pm.enu2aer(te,tn,tu), (e2az,e2el,e2rn),
                    rtol=0.01,
                    err_msg='enu2aer: ' + str(pm.enu2aer(te,tn,tu)))

    assert_allclose(pm.ned2aer(tn,te,-tu), (e2az,e2el,e2rn),
                    rtol=0.01,
                    err_msg='enu2aer: ' + str(pm.enu2aer(te,tn,tu)))

    assert_allclose(pm.enu2geodetic(te,tn,tu,tlat,tlon,talt),(lat2,lon2,alt2),
                    rtol=0.01,
                    err_msg='enu2geodetic: '+str(pm.enu2geodetic(te,tn,tu,tlat,
                                                              tlon,talt)))

    assert_allclose(pm.ned2geodetic(tn,te,-tu,tlat,tlon,talt),(lat2,lon2,alt2),
                    rtol=0.01,
                    err_msg='enu2geodetic: '+str(pm.enu2geodetic(te,tn,tu,tlat,
                                                              tlon,talt)))

    assert_allclose(pm.enu2ecef(te,tn,tu,tlat,tlon,talt),(e2x, e2y, e2z),
                    rtol=0.01,
                    err_msg='enu2ecef: '+
                    str(pm.enu2ecef(te,tn,tu,tlat,tlon,talt)))

    assert_allclose(pm.ned2ecef(tn,te,-tu,tlat,tlon,talt),(e2x,e2y,e2z),
                    rtol=0.01,
                    err_msg='ned2ecef: '+
                    str(pm.ned2ecef(tn,te,-tu,tlat,tlon,talt)))
#%%
def test_eci():
    tlla = (42,-82,200)
    teci = (-3.977913815668146e6,-2.582332196263046e6,4.250818828152067e6)
    t = datetime(2013,1,15,12,0,5,tzinfo=UTC)
    lla = asarray(pm.eci2geodetic(teci,t)).squeeze()
    assert_allclose(lla,tlla,rtol=0.2)

    assert_allclose(pm.eci2ecef(teci,t).squeeze(),
                [649012.04640917,-4697980.55129606,4250818.82815207])

    assert_allclose(pm.ecef2eci([649012.04640917,-4697980.55129606,4250818.82815207],t).squeeze(),
                    teci)

    assert_allclose(asarray(pm.eci2aer(teci,42,-100,0,t)).squeeze(),
                    [83.73050,-6.614478,1.473510e6])

def test_vdist():
    dist_m = vdist(10,20,10.021372672660874,20.016847098929979)[0]
    assert_almost_equal(dist_m,2999.9999971763464)
#%%
if __name__ == '__main__':
    run_module_suite()
