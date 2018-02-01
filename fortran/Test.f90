program Test

use iso_fortran_env, only: wp=>real64
use maptran

implicit none

real(wp), parameter :: lat = 42, lon= -82, alt = 200, &
                       az = 33, el=70 ,srange = 1e3, &
                       x0 = 660.675e3, y0 = -4700.949e3, z0 = 4245.738e3, &
                       xl = 660.930e3, yl = -4701.424e3, zl = 4246.579e3, & ! aer2ecef
                       er = 186.277521, nr = 286.84222, ur = 939.69262 ! aer2enu

real(wp) :: lat2, lon2, alt2,x1,y1,z1,x2,y2,z2,x3,y3,z3,az3,el3,rng3,e1,n1,u1,e3,n3,u3

type(wgs84Ellipsoid) :: ell

ell = wgs84Ellipsoid()


call ecef2geodetic(x0,y0,z0,lat2,lon2,alt2)
call assert_isclose(lat2,lat)
call assert_isclose(lon2,lon)
call assert_isclose(alt2,alt)

call geodetic2ecef(lat,lon,alt,x1,y1,z1)
call assert_isclose(x1,x0)
call assert_isclose(y1,y0)
call assert_isclose(z1,z0)

call aer2enu(az, el, srange, e1, n1, u1)
call assert_isclose(e1, er)

call aer2ecef(az,el,srange,lat,lon,alt,x2,y2,z2)
call assert_isclose(x2, xl)
call assert_isclose(y2, yl)
call assert_isclose(z2, zl)

call enu2ecef(e1,n1,u1,lat,lon,alt, x3, y3, z3)
call assert_isclose(x3, x2)

call ecef2enu(x3,y3,z3,lat,lon,alt, e3,n3,u3)
call assert_isclose(e3,e1)

call ecef2aer(x2,y2,z2, lat,lon,alt, az3, el3, rng3)
call assert_isclose(el3, el)
call assert_isclose(az3, az)


end program
