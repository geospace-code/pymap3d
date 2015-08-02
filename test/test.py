#!/usr/bin/env python
"""
runs tests
"""
from __future__ import absolute_import, division
from numpy.testing import assert_allclose, assert_almost_equal
#
from pymap3d.coordconv3d import *
from pymap3d.azel2radec import azel2radec
from pymap3d.haversine import angledist,angular_separation
from pymap3d.vreckon import vreckon


def test_azel2radec():
    ra,dec = azel2radec(180.1, 80, 65, -148, '2014-04-06T08:00:00Z')
    assert_allclose(ra,166.5032081149338,rtol=1e-2)
    assert_allclose(dec,55.000011165405752,rtol=1e-2)

def test_haversine():
    assert_almost_equal(angledist(35,23,84,20),45.482789587392013)
        #%% compare with astropy
    dist_astropy = angular_separation(radians(35),radians(23),
                                      radians(84),radians(20))
    assert_almost_equal(45.482789587392013, dist_astropy)

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

#enu2ccef
e2x, e2y, e2z = (6.679456743004259e+05, -4.693230928738789e+06, 4.254723326333052e+06)

e2az, e2el, e2rn = (36.664402767128749, 4.477194667550686, 1.389837889201037e+04)


def test_geodetic():
    assert_allclose(ecef2geodetic(tx,ty,tz),(ec2la,ec2lo,ec2a),
                rtol=0.01,
                err_msg='ecef2geodetic: ' + str(ecef2geodetic(tx,ty,tz))
                )

    assert_allclose(geodetic2aer(lat2,lon2,alt2,tlat,tlon,talt),
                    (g2az,g2el,g2rn),
                    rtol=0.05,
                    err_msg= 'geodetic2aer: {}'.format(
                    geodetic2aer(lat2,lon2,alt2,tlat,tlon,talt)))

    assert_allclose(geodetic2ecef(tlat,tlon,talt),(g2x,g2y,g2z),
                    rtol=0.01,
                    err_msg='geodetic2ecef: {}'.format(
                    geodetic2ecef(tlat,tlon,talt)))

def test_ecefenu():
    assert_allclose(aer2ecef(taz,tel,tsrange,tlat,tlon,talt), (a2x,a2y,a2z),
                    rtol=0.01,
                    err_msg='aer2ecef: {}'.format(
                    aer2ecef(taz,tel,tsrange,tlat,tlon,talt)))

    assert_allclose(aer2enu(taz,tel,tsrange),(a2e,a2n,a2u),
                    rtol=0.01,
                    err_msg='aer2enu: ' + str(aer2enu(taz,tel,tsrange)))

    assert_allclose(aer2ned(taz,tel,tsrange),(a2n,a2e,-a2u),
                    rtol=0.01,
                    err_msg='aer2ned: ' + str(aer2ned(taz,tel,tsrange)))


    assert_allclose(ecef2enu(tx,ty,tz, tlat, tlon, talt),(e2e,e2n,e2u),
                    rtol=0.01,
                    err_msg='ecef2enu: {}'.format(ecef2enu(tx,ty,tz, tlat,
                    tlon, talt)))

    assert_allclose(ecef2ned(tx,ty,tz, tlat, tlon, talt),(e2n,e2e,-e2u),
                    rtol=0.01,
                    err_msg='ecef2ned: {}'.format(ecef2enu(tx,ty,tz, tlat,
                    tlon, talt)))

    assert_allclose(aer2geodetic(taz,tel,tsrange,tlat,tlon,talt),(a2la,a2lo,a2a),
                rtol=0.01,err_msg='aer2geodetic {}'.format(
                  aer2geodetic(taz,tel,tsrange,tlat,tlon,talt)))

    assert_allclose(ecef2aer(tx, ty, tz, tlat, tlon,talt),
                    (ec2az,ec2el,ec2rn),
                    rtol=0.01,
                    err_msg='ecef2aer {}'.format(ecef2aer(a2x, a2y, a2z, tlat,
                    tlon, talt)))

    assert_allclose(enu2aer(te,tn,tu), (e2az,e2el,e2rn),
                    rtol=0.01,
                    err_msg='enu2aer: ' + str(enu2aer(te,tn,tu)))

    assert_allclose(ned2aer(tn,te,-tu), (e2az,e2el,e2rn),
                    rtol=0.01,
                    err_msg='enu2aer: ' + str(enu2aer(te,tn,tu)))

    assert_allclose(enu2geodetic(te,tn,tu,tlat,tlon,talt),(lat2,lon2,alt2),
                    rtol=0.01,
                    err_msg='enu2geodetic: '+str(enu2geodetic(te,tn,tu,tlat,
                                                              tlon,talt)))

    assert_allclose(ned2geodetic(tn,te,-tu,tlat,tlon,talt),(lat2,lon2,alt2),
                    rtol=0.01,
                    err_msg='enu2geodetic: '+str(enu2geodetic(te,tn,tu,tlat,
                                                              tlon,talt)))

    assert_allclose(enu2ecef(te,tn,tu,tlat,tlon,talt),(e2x, e2y, e2z),
                    rtol=0.01,
                    err_msg='enu2ecef: '+
                    str(enu2ecef(te,tn,tu,tlat,tlon,talt)))

    assert_allclose(ned2ecef(tn,te,-tu,tlat,tlon,talt),(e2x,e2y,e2z),
                    rtol=0.01,
                    err_msg='ned2ecef: '+
                    str(ned2ecef(tn,te,-tu,tlat,tlon,talt)))

    assert_allclose(eci2ecef((10e6,20e6,30e6),.230).squeeze(),
                (1.429621442075752e7,1.719355266475562e7,3e7),
                rtol=0.01)

    assert_allclose(ecef2eci((1.429621442075752e7,1.719355266475562e7,3e7),.230).squeeze(),
                (10e6,20e6,30e6),
                rtol=0.01)

if __name__ == '__main__':
    test_azel2radec()
    test_haversine()
    test_vreckon()

    test_geodetic()
    test_ecefenu()