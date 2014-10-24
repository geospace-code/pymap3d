#!/usr/bin/env python3
# Michael Hirsch ported and adaptation from
# GNU Octave Mapping Toolbox by
# Copyright (c) 2013, Sandeep V. Mahanthi
# Copyright (c) 2013, Felipe G. Nievinski
from __future__ import division
from numpy import sin,cos,tan,sqrt,radians,arctan2,hypot,degrees,mod,atleast_2d,empty_like,array


class EarthEllipsoid:
    def __init__(self):
        self.a = 6378137.0  # semi-major axis [m]
        self.f = 1.0 / 298.2572235630  # flattening
        self.b = self.a * (1 - self.f)  # semi-minor axis
        self.e = sqrt ( (self.a**2 - self.b**2) / (self.a**2)) # first eccentricity

# Parameters needed to evaluate the geopotential.
# I got them from the US National Geospatial-Inteligence Agency's
# website, "NGA/NASA EGM96, N=M=360 Earth Gravity Model"
# (NGA > Products and Services >Geospatial Sciences Division
# > Physical Geodesy Home)
# <http://earth-info.nga.mil/GandG/wgsegm/egm96.html>
#
# Earth's Gravitational Constant w/ atmosphere:
        self.GM = 0.3986004418e15  # m^3/s^2

        # Earth's angular speed:
        self.omega = 7292115.0e-11  # rad/s

# Dynamical form factor of the Earth:
# The four defining paramenters of WGS-84 are
# GM, omega, a, and 1/f. We need J2 instead of f:
       # self.J2 = get_J2 (self.a, self.b, self.omega, self.GM)

def aer2ecef(az,el,srange,lat0,lon0,alt0,ell=EarthEllipsoid()):
    # Origin of the local system in geocentric coordinates.
    x0,y0,z0 = geodetic2ecef(lat0,lon0,alt0,ell)
    # Convert Local Spherical AER to ENU
    e1, n1, u1 = aer2enu(az, el, srange)
    # Rotating ENU to ECEF
    dx, dy, dz = enu2ecef_int(e1, n1, u1, lat0, lon0)
    # Origin + offset from origin equals position in ECEF
    return x0 + dx, y0 + dy, z0 + dz

def aer2enu(az,el,srange):
    #Calculation of AER2ENU
    u1 = srange * sin(radians(el))
    r   = srange * cos(radians(el))
    e1  = r * sin(radians(az))
    n1 = r * cos(radians(az))

    return e1,n1,u1

def aer2geodetic(az,el,srange,lat0,lon0,alt0):
#angles in DEGREES
#range in METERS
    x,y,z = aer2ecef(az,el,srange,lat0,lon0,alt0)
    return ecef2geodetic(x,y,z)

def aer2ned(az,elev,slantRange):
    xNorth,yEast,zUp = aer2enu(az,elev,slantRange)
    zDown = -zUp
    return xNorth,yEast,zDown

def ecef2aer(x, y, z, lat0, lon0, h0, ell=EarthEllipsoid()):
    xEast, yNorth, zUp = ecef2enu(x, y, z, lat0, lon0, h0, ell)
    return enu2aer(xEast, yNorth, zUp)

def ecef2enu(x, y, z, lat0, lon0, h0, ell=EarthEllipsoid()):
    x0,y0,z0 = geodetic2ecef(lat0, lon0, h0,ell)
    xEast, yNorth, zUp = ecef2enu_int(x - x0, y - y0, z - z0, lat0, lon0)
    return xEast,yNorth,zUp

def ecef2enu_int(u, v, w, lat0, lon0):
    t     =  cos(radians(lon0)) * u + sin(radians(lon0)) * v
    East  = -sin(radians(lon0)) * u + cos(radians(lon0)) * v
    Up    =  cos(radians(lat0)) * t + sin(radians(lat0)) * w
    North = -sin(radians(lat0)) * t + cos(radians(lat0)) * w
    return East, North, Up

def ecef2geodetic(x,y=None,z=None,ell=EarthEllipsoid()):
    if y is None:
        x,y,z = depack(x)
    #Algorithm is based on
    #http://www.astro.uni.torun.pl/~kb/Papers/geod/Geod-BG.htm
    #This algorithm provides a converging solution to the latitude equation
    #in terms of the parametric or reduced latitude form (v)
    #This algorithm provides a uniform solution over all latitudes as it does
    #not involve division by cos(phi) or sin(phi)
    a = ell.a
    b = ell.b
    r = hypot(x,y)
# Constant required for Latitude equation
    rho = arctan2(b*z,a*r)
#Constant required for latitude equation
    c = (a**2 - b**2) / sqrt((a*r)**2 + (b*z)**2)
# Starter for the Newtons Iteration Method
    vnew = arctan2(a*z,b*r)
# Initializing the parametric latitude
    v = 0.
    count = 0
    while (v != vnew).any() and count<5:
        v = vnew
# Derivative of latitude equation
        w = 2*(cos(v - rho) - c*cos(2.*v))
# Newtons Method for computing iterations
        vnew = v - ((2*sin(v-rho)-c*sin(2.*v))/w)
 #       print(count)
        count+=1

    #Computing latitude from the root of the latitude equation
    lat = arctan2((a*tan(vnew)),b)
    #by inspection
    lon = arctan2(y,x)

    alt = ((r-(a*cos(vnew)))*cos(lat)) + ((z-(b*sin(vnew)))*sin(lat))

    #assert all(-90<=lat) and all(lat<=90)
    #assert all(-180<=lon) and all(lon<=180)

    return degrees(lat),degrees(lon),alt

def ecef2ned(x, y, z, lat0, lon0, h0, ell=EarthEllipsoid()):
    yEast, xNorth, zUp = ecef2enu(x, y, z, lat0, lon0, h0, ell)
    zDown = -zUp
    return xNorth, yEast, zDown

def ecef2ned_int(x, y, z, lat0, lon0):
    xEast, yNorth, zUp = ecef2enu_int(x, y, z, lat0, lon0)
    zDown = -zUp
    return xEast, yNorth, zDown

def enu2aer(xEast, yNorth, zUp):
    r = hypot(xEast, yNorth)
    slantRange = hypot(r, zUp)
    elev = degrees(arctan2(zUp,r))
    az = degrees(mod(arctan2(xEast,yNorth), 2*arctan2(0,-1)))
    return az, elev, slantRange

def enu2ecef(e1,n1,u1,lat0,lon0,alt0,ell):
    x0, y0, z0 = geodetic2ecef(lat0, lon0, alt0, ell)
    dx, dy, dz = enu2ecef_int(e1, n1, u1, lat0, lon0)
    return x0 + dx, y0 + dy, z0 + dz

def enu2ecef_int(es,nr,up,lat0,lon0):
    t = cos(radians(lat0)) * up - sin(radians(lat0)) * nr
    w = sin(radians(lat0)) * up + cos(radians(lat0)) * nr

    u = cos(radians(lon0)) * t - sin(radians(lon0)) * es
    v = sin(radians(lon0)) * t + cos(radians(lon0)) * es
    return u,v,w

def enu2geodetic(East, North, Up, lat0, lon0, h0, ell):
    x, y, z = enu2ecef(East, North, Up, lat0, lon0, h0, ell)
    return ecef2geodetic(x, y, z, ell)
#%% These are ported -- have not carefully verified yet.
def eci2ecef(eci,lst):
    eci = atleast_2d(eci)
    N,trip = eci.shape
    if eci.ndim > 2 or trip != 3:
        raise RuntimeError('eci triplets must be shape (N,3)')
    '''
    ported from:
    https://github.com/dinkelk/astrodynamics/blob/master/rot3.m
    '''
    ecef = empty_like(eci)
    for i in range(N):
        ecef[i,:] = rottrip(lst[i]).dot(eci[i,:])
    return ecef

def ecef2eci(ecef,lst):
    ecef = atleast_2d(ecef)
    N,trip = ecef.shape
    if ecef.ndim > 2 or trip != 3:
        raise RuntimeError('ecef triplets must be shape (N,3)')
    '''
    ported from:
    https://github.com/dinkelk/astrodynamics/blob/master/rot3.m
    '''
    eci = empty_like(ecef)
    for i in range(N):
        eci[i,:] = rottrip(lst[i]).T.dot(ecef[i,:]) #this one is transposed
    return eci

def rottrip(ang):
    ang = ang.squeeze()
    if ang.size>1:
        raise RuntimeError('this function is for one angle at a time')
    '''
    ported from:
    https://github.com/dinkelk/astrodynamics/blob/master/rot3.m
    '''
    return array([[cos(ang),  sin(ang), 0],
                 [-sin(ang), cos(ang), 0],
                 [0,         0,        1]])
#%%
def geodetic2aer(lat, lon, h, lat0, lon0, h0, ell=EarthEllipsoid()):
    East, North, Up = geodetic2enu(lat, lon, h, lat0, lon0, h0, ell)
    return enu2aer(East, North, Up)


def geodetic2ecef(lat,lon,alt,ell=EarthEllipsoid()):
# Auxiliary quantities
# radius of curvature of the prime vertical section
    N = compute_prime_vertical_radius (lat, ell)
    # Some shortnames for variables used often.
    a = ell.a;  b = ell.b

    # Compute cartesian (geocentric) coordinates given  (curvilinear) geodetic coordinates.
    x = (N + alt) * cos(radians(lat))  * cos(radians(lon))
    y = (N + alt) * cos(radians(lat))  * sin(radians(lon))
    z = (N * (b/a)**2 + alt) * sin(radians(lat))
    return x,y,z

def geodetic2enu(lat, lon, h, lat0, lon0, h0, ell=EarthEllipsoid()):
    x1,y1,z1 = geodetic2ecef(lat,lon,h,ell)
    x2,y2,z2 = geodetic2ecef(lat0,lon0,h0,ell)
    dx = x1-x2
    dy = y1-y2
    dz = z1-z2
    return ecef2enu_int(dx, dy, dz, lat0, lon0)

def geodetic2ned(lat, lon, h, lat0, lon0, h0, ell=EarthEllipsoid()):
    yEast, xNorth, zUp = geodetic2enu(lat, lon, h, lat0, lon0, h0, ell)
    return xNorth, yEast, -zUp

def ned2aer(xNorth, yEast, zDown):
    return enu2aer(yEast, xNorth, -zDown)

def ned2ecef(xNorth, yEast, zDown, lat0, lon0, h0, ell=EarthEllipsoid()):
    return enu2ecef(yEast, xNorth, -zDown, lat0, lon0, h0, ell)

def ned2ecef_int(uNorth, vEast, wDown, lat0, lon0):
    return enu2ecef_int(vEast, uNorth, -wDown, lat0, lon0)

def ned2geodetic(xNorth, yEast, zDown, lat0, lon0, h0, ell=EarthEllipsoid()):
    x, y, z = enu2ecef(yEast, xNorth, -zDown, lat0, lon0, h0, ell)
    return ecef2geodetic(x, y, z,ell)


def compute_prime_vertical_radius (lat, ell):
    N = get_radius_normal (lat, ell)
    return N

def get_radius_normal(lat,ell):
    a = ell.a;  b = ell.b;
    lat2 = radians(lat)
    return a**2 / sqrt( a**2 * (cos(lat2))**2 + b**2 * (sin(lat2))**2 )

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
    #test suite
    lat,lon,alt = 42., -82., 200.
    az,el,srange = 45., 80., 1000.
    x, y, z = 23893., 2193., 2389239.
#test results
    et,nt,ut = aer2enu(az,el,srange)
    xt,yt,zt = aer2ecef(az,el,srange,lat,lon,alt)
    latt,lont,altt = aer2geodetic(az,el,srange,lat,lon,alt)
    azt,elt,rngt = ecef2aer(x, y, z, lat, lon, alt)

    print('aer2enu ' + str((et,nt,ut)))
    print('aer2ecef ' + str((xt,yt,zt)))
    print('aer2geodetic ' + str((latt,lont,altt)))
    print('ecef2aer ' + str((azt,elt,rngt)))