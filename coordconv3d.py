"""Michael Hirsch ported and adaptation from
 GNU Octave Mapping Toolbox by
  Copyright (c) 2013, Sandeep V. Mahanthi
 Copyright (c) 2013, Felipe G. Nievinski

 Input/output: units are METERS and DEGREES. boolean deg=True means degrees
"""
from __future__ import division
from numpy import (sin,cos,tan,sqrt,radians,arctan2,hypot,degrees,mod,
                   atleast_2d,atleast_1d,empty_like,array)

class EarthEllipsoid:
    def __init__(self):
        self.a = 6378137.0  # semi-major axis [m]
        self.f = 1.0 / 298.2572235630  # flattening
        self.b = self.a * (1 - self.f)  # semi-minor axis

def aer2ecef(az,el,srange,lat0,lon0,alt0,ell=EarthEllipsoid(),deg=True):
    # Origin of the local system in geocentric coordinates.
    x0,y0,z0 = geodetic2ecef(lat0,lon0,alt0,ell,deg=deg)
    # Convert Local Spherical AER to ENU
    e1, n1, u1 = aer2enu(az, el, srange, deg=deg)
    # Rotating ENU to ECEF
    dx, dy, dz = enu2ecef_int(e1, n1, u1, lat0, lon0,deg=deg)
    # Origin + offset from origin equals position in ECEF
    return x0 + dx, y0 + dy, z0 + dz

def aer2enu(az,el,srange,deg=True):
    if deg:
        el = radians(el)
        az = radians(az)

    r   = srange * cos(el)
    return r*sin(az), r*cos(az), srange*sin(el)

def aer2geodetic(az,el,srange,lat0,lon0,alt0,deg=True):
    x,y,z = aer2ecef(az,el,srange,lat0,lon0,alt0,deg=deg)
    return ecef2geodetic(x,y,z,deg=deg)

def aer2ned(az,elev,slantRange,deg=True):
    xNorth,yEast,zUp = aer2enu(az,elev,slantRange,deg=deg)
    return xNorth,yEast,-zUp

def ecef2aer(x, y, z, lat0, lon0, h0, ell=EarthEllipsoid(),deg=True):
    xEast, yNorth, zUp = ecef2enu(x, y, z, lat0, lon0, h0, ell,deg=deg)
    return enu2aer(xEast, yNorth, zUp,deg=deg)

def ecef2enu(x, y, z, lat0, lon0, h0, ell=EarthEllipsoid(),deg=True):
    x0,y0,z0 = geodetic2ecef(lat0, lon0, h0,ell,deg=deg)
    return ecef2enu_int(x - x0, y - y0, z - z0, lat0, lon0,deg=deg)

def ecef2enu_int(u, v, w, lat0, lon0,deg=True):
    if deg:
        lat0 = radians(lat0)
        lon0 = radians(lon0)
    t     =  cos(lon0) * u + sin(lon0) * v
    East  = -sin(lon0) * u + cos(lon0) * v
    Up    =  cos(lat0) * t + sin(lat0) * w
    North = -sin(lat0) * t + cos(lat0) * w
    return East, North, Up

def ecef2geodetic(x,y=None,z=None,ell=EarthEllipsoid(),deg=True):
    if y is None:
        x,y,z = depack(x)
    """Algorithm is based on
    http://www.astro.uni.torun.pl/~kb/Papers/geod/Geod-BG.htm
    This algorithm provides a converging solution to the latitude equation
    in terms of the parametric or reduced latitude form (v)
    This algorithm provides a uniform solution over all latitudes as it does
    not involve division by cos(phi) or sin(phi)
    """
    ea = ell.a
    eb = ell.b
    rad = hypot(x,y)
# Constant required for Latitude equation
    rho = arctan2(eb*z,ea*rad)
#Constant required for latitude equation
    c = (ea**2 - eb**2) / hypot(ea*rad,eb*z)
# Starter for the Newtons Iteration Method
    vnew = arctan2(ea*z,eb*rad)
# Initializing the parametric latitude
    v = 0
    count = 0
    while (v != vnew).any() and count<5:
        v = vnew.copy()
# Newtons Method for computing iterations
        vnew = v - ((2*sin(v-rho)-c*sin(2*v)) /
                    (2*(cos(v - rho) - c*cos(2*v))) )
 #       print(count)
        count+=1

    #Computing latitude from the root of the latitude equation
    lat = arctan2(ea*tan(vnew),eb)
    #by inspection
    lon = arctan2(y,x)

    alt = ((rad-ea*cos(vnew))*cos(lat)) + ((z-eb*sin(vnew))*sin(lat))

    if deg:
        return degrees(lat),degrees(lon),alt
    else:
        return lat, lon, alt #radians

# this is from PySatel and give virtually an identical result to EIGHT decimal places!
#def cbrt(x):
#	if x >= 0:
#		return pow(x, 1.0/3.0)
#	else:
#		return -pow(abs(x), 1.0/3.0)
#
#def ecef2geodetic(x, y, z, ell=EarthEllipsoid(),deg=True):
#    a = ell.a; b = ell.b
#    esq = 6.69437999014*0.001
#    e1sq = 6.73949674228*0.001
#    r = hypot(x,y)
#    Esq = a**2 - b**2
#    F = 54 * b**2 * z**2
#    G = r**2 + (1 - esq)* z**2 - esq*Esq
#    C = (esq**2 *F* r**2)/(pow(G, 3))
#    S = cbrt(1 + C + sqrt(C**2 + 2*C))
#    P = F/(3* pow((S + 1/S + 1), 2)*G**2)
#    Q = sqrt(1 + 2* esq**2 *P)
#    r_0 =  -(P*esq*r)/(1 + Q) + sqrt(0.5* a**2 *(1 + 1.0/Q) - \
#           P*(1 - esq)*z**2/(Q*(1 + Q)) - 0.5*P* r**2)
#    U = sqrt(pow((r - esq*r_0), 2) + z**2)
#    V = sqrt(pow((r - esq*r_0), 2) + (1 - esq)* z**2)
#    Z_0 = b**2 *z/(a*V)
#    alt = U*(1 - b**2/(a*V))
#    lat = arctan((z + e1sq*Z_0)/r)
#    lon = arctan2(y, x)
#
#    if deg:
#        return degrees(lat),degrees(lon),alt
#    else:
#        return lat, lon, alt #radians

def ecef2ned(x, y, z, lat0, lon0, h0, ell=EarthEllipsoid(),deg=True):
    yEast, xNorth, zUp = ecef2enu(x, y, z, lat0, lon0, h0, ell,deg=deg)
    return xNorth, yEast, -zUp

def ecef2ned_int(x, y, z, lat0, lon0,deg=True):
    xEast, yNorth, zUp = ecef2enu_int(x, y, z, lat0, lon0,deg=deg)
    return xEast, yNorth, -zUp

def enu2aer(e, n, u,deg=True):
    r = hypot(e, n)
    slantRange = hypot(r, u)
    elev = arctan2(u,r)
    az = mod(arctan2(e,n), 2*arctan2(0,-1))
    if deg:
        return degrees(az), degrees(elev), slantRange
    else:
        return az, elev, slantRange #radians

def enu2ecef(e1,n1,u1,lat0,lon0,alt0,ell=EarthEllipsoid(),deg=True):
    x0, y0, z0 = geodetic2ecef(lat0, lon0, alt0, ell,deg=deg)
    dx, dy, dz = enu2ecef_int(e1, n1, u1, lat0, lon0,deg=deg)
    return x0 + dx, y0 + dy, z0 + dz

def enu2ecef_int(east,north,up,lat0,lon0,deg=True):
    if deg:
        lat0 = radians(lat0)
        lon0 = radians(lon0)
    t = cos(lat0) * up - sin(lat0) * north
    w = sin(lat0) * up + cos(lat0) * north

    u = cos(lon0) * t - sin(lon0) * east
    v = sin(lon0) * t + cos(lon0) * east
    return u,v,w

def enu2geodetic(e, n, u, lat0, lon0, h0, ell=EarthEllipsoid(),deg=True):
    x, y, z = enu2ecef(e, n, u, lat0, lon0, h0, ell,deg=deg)
    return ecef2geodetic(x, y, z, ell,deg=deg)
#====================================================================
#%%
    """inputs:

    ece/ecef: a Nx3 vector of x,y,z triplets in the eci or ecef system [meters]
    lst: length N vector of sidereal time angle [radians]. The function datetime2hourangle.py in
    https://github.com/scienceopen/astrometry can provide this for you.
    """
def eci2ecef(eci,lst):
    lst = atleast_1d(lst)
    eci = atleast_2d(eci)
    N,trip = eci.shape
    if eci.ndim > 2 or trip != 3:
        exit('eci2ecef: eci triplets must be shape (N,3)')
    """ported from:
    https://github.com/dinkelk/astrodynamics/blob/master/rot3.m
    """
    ecef = empty_like(eci)
    for i in range(N):
        ecef[i,:] = rottrip(lst[i]).dot(eci[i,:])
    return ecef

def ecef2eci(ecef,lst):
    lst = atleast_1d(lst)
    ecef = atleast_2d(ecef)
    N,trip = ecef.shape
    if ecef.ndim > 2 or trip != 3:
        exit('ecef2eci: ecef triplets must be shape (N,3)')
    """ported from:
    https://github.com/dinkelk/astrodynamics/blob/master/rot3.m
    """
    eci = empty_like(ecef)
    for i in range(N):
        eci[i,:] = rottrip(lst[i]).T.dot(ecef[i,:]) #this one is transposed
    return eci

def rottrip(ang):
    ang = ang.squeeze()
    if ang.size>1:
        exit('rottrip: only one angle allowed at a time')
    """ported from:
    https://github.com/dinkelk/astrodynamics/blob/master/rot3.m
    """
    return array([[cos(ang),  sin(ang), 0],
                 [-sin(ang), cos(ang), 0],
                 [0,         0,        1]])
#==========================================================
#%%
def geodetic2aer(lat, lon, h, lat0, lon0, h0, ell=EarthEllipsoid(),deg=True):
    e, n, u = geodetic2enu(lat, lon, h, lat0, lon0, h0, ell,deg=deg)
    return enu2aer(e, n, u,deg=deg)

def geodetic2ecef(lat,lon,alt,ell=EarthEllipsoid(),deg=True):
    if deg:
        lat = radians(lat)
        lon = radians(lon)
    # radius of curvature of the prime vertical section
    N = get_radius_normal(lat, ell)
    # Compute cartesian (geocentric) coordinates given  (curvilinear) geodetic coordinates.
    x = (N + alt) * cos(lat)  * cos(lon)
    y = (N + alt) * cos(lat)  * sin(lon)
    z = (N * (ell.b/ell.a)**2 + alt) * sin(lat)
    return x,y,z

def geodetic2enu(lat, lon, h, lat0, lon0, h0, ell=EarthEllipsoid(),deg=True):
    x1,y1,z1 = geodetic2ecef(lat,lon,h,ell,deg=deg)
    x2,y2,z2 = geodetic2ecef(lat0,lon0,h0,ell,deg=deg)
    dx = x1-x2
    dy = y1-y2
    dz = z1-z2
    return ecef2enu_int(dx, dy, dz, lat0, lon0,deg=deg)

def geodetic2ned(lat, lon, h, lat0, lon0, h0, ell=EarthEllipsoid(),deg=True):
    yEast, xNorth, zUp = geodetic2enu(lat, lon, h, lat0, lon0, h0, ell,deg=deg)
    return xNorth, yEast, -zUp

def ned2aer(xNorth, yEast, zDown,deg=True):
    return enu2aer(yEast, xNorth, -zDown,deg=deg)

def ned2ecef(xNorth, yEast, zDown, lat0, lon0, h0, ell=EarthEllipsoid(),deg=True):
    return enu2ecef(yEast, xNorth, -zDown, lat0, lon0, h0, ell,deg=deg)

def ned2ecef_int(uNorth, vEast, wDown, lat0, lon0,deg=True):
    return enu2ecef_int(vEast, uNorth, -wDown, lat0, lon0,deg=deg)

def ned2geodetic(xNorth, yEast, zDown, lat0, lon0, h0, ell=EarthEllipsoid(),deg=True):
    x, y, z = enu2ecef(yEast, xNorth, -zDown, lat0, lon0, h0, ell,deg=deg)
    return ecef2geodetic(x, y, z,ell,deg=deg)


def get_radius_normal(lat_radians,ell):
    a = ell.a;  b = ell.b
    return a**2 / sqrt( a**2 * (cos(lat_radians))**2 + b**2 * (sin(lat_radians))**2 )

def depack(x0):
    if x0.ndim>2:
        raise RuntimeError('I expect Nx3 or 3XN triplets')
    m,n = x0.shape
    if m == 3: # 3xN triplets
        x = x0[0,:]
        y = x0[1,:]
        z = x0[2,:]
    elif n==3: # Nx3 triplets
        x = x0[:,0]
        y = x0[:,1]
        z = x0[:,2]
    else:
        raise RuntimeError('I expect an Nx3 or 3xN input of x,y,z')
    return x,y,z

if __name__ == '__main__':
    #selftest
    from numpy.testing import assert_allclose
    #test suite
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

    g2x, g2y, g2z = 660675.251825, -4700948.683162, 4245737.662222 #geodeteic2ecef
    e2e, e2n, e2u = 8272.476048, 11112.773942, 84.941624 #ecef2enu

    ec2la, ec2lo, ec2a = 42.100000, -81.900000, 300.000000 #ecef2geodetic

    e2az, e2el, e2rn = 36.664402767128749, 4.477194667550686, 1.389837889201037e+04
    e2x, e2y, e2z = 2.198984328830889e+06, -1.084794996374469e+07, 3.605050273624581e+06
#test results
    assert_allclose(ecef2geodetic(tx,ty,tz),(ec2la,ec2lo,ec2a),
                    rtol=0.01,
                    err_msg='ecef2geodetic: ' + str(ecef2geodetic(tx,ty,tz)) )

    assert_allclose(geodetic2aer(lat2,lon2,alt2,tlat,tlon,talt), (g2az,g2el,g2rn),
                    rtol=0.05,
                    err_msg= 'geodetic2aer: ' + str(geodetic2aer(lat2,lon2,alt2,tlat,tlon,talt)))

    assert_allclose(geodetic2ecef(tlat,tlon,talt),(g2x,g2y,g2z),
                    rtol=0.01,
                    err_msg='geodetic2ecef: ' + str(geodetic2ecef(tlat,tlon,talt)))

    assert_allclose(aer2ecef(taz,tel,tsrange,tlat,tlon,talt), (a2x,a2y,a2z),
                             rtol=0.01,
                             err_msg='aer2ecef: ' + str(aer2ecef(taz,tel,tsrange,tlat,tlon,talt)))

    assert_allclose(aer2enu(taz,tel,tsrange),(a2e,a2n,a2u),
                    rtol=0.01,
                    err_msg='aer2enu: ' + str(aer2enu(taz,tel,tsrange)))

    assert_allclose(ecef2enu(tx,ty,tz, tlat, tlon, talt),(e2e,e2n,e2u),
                    rtol=0.01,
                    err_msg='ecef2enu: ' + str(ecef2enu(tx,ty,tz, tlat, tlon, talt)))

    assert_allclose(aer2geodetic(taz,tel,tsrange,tlat,tlon,talt),(a2la,a2lo,a2a),
                    rtol=0.01,err_msg='aer2geodetic' + str(aer2geodetic(taz,tel,tsrange,tlat,tlon,talt)))

    assert_allclose(ecef2aer(tx, ty, tz, tlat, tlon,talt), (ec2az,ec2el,ec2rn),
                    rtol=0.01,
                    err_msg='ecef2aer' + str(ecef2aer(a2x, a2y, a2z, tlat, tlon, talt)))

    assert_allclose(enu2aer(te,tn,tu), (e2az,e2el,e2rn),
                    rtol=0.01,
                    err_msg='enu2aer: ' + str(enu2aer(te,tn,tu)))

    assert_allclose(enu2geodetic(te,tn,tu,tlat,tlon,talt),(lat2,lon2,alt2),
                    rtol=0.01,
                    err_msg='enu2geodetic: ' + str(enu2geodetic(te,tn,tu,tlat,tlon,talt)))

    assert_allclose(enu2ecef(tx,ty,tz,tlat,tlon,talt),(e2x,e2y,e2z),
                    rtol=0.01,
                    err_msg='enu2ecef: '+ str(enu2ecef(tx,ty,tz,tlat,tlon,talt)))
    assert_allclose(eci2ecef((10e6,20e6,30e6),.230).squeeze(),
                    (1.429621442075752e7,1.719355266475562e7,3e7),
                    rtol=0.01)
    assert_allclose(ecef2eci((1.429621442075752e7,1.719355266475562e7,3e7),.230).squeeze(),
                    (10e6,20e6,30e6),
                    rtol=0.01)

    exit(0)