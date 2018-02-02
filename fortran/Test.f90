program Test

use iso_fortran_env, only: wp=>real64
use maptran

implicit none

real(wp), parameter :: lat = 42, lon= -82, alt = 200, &
                       az = 33, el=70, rng= 1e3, &
                       x0 = 660.675e3, y0 = -4700.949e3, z0 = 4245.738e3, &
                       xl = 660.930e3, yl = -4701.424e3, zl = 4246.579e3, & ! aer2ecef
                       er = 186.277521, nr = 286.84222, ur = 939.69262, & ! aer2enu
                       lat1 = 42.0026, lon1 = -81.9978, alt1 = 1.1397e3 ! aer2geodetic

real(wp) :: lat2, lon2, alt2,lat3,lon3,alt3,lat4,lon4,alt4,&
            x1,y1,z1,x2,y2,z2,x3,y3,z3,&
            az2,el2,rng2,az3,el3,rng3,az4,el4,rng4,&
            e1,n1,u1,e2,n2,u2,e3,n3,u3

type(referenceEllipsoid), parameter :: spheroid = wgs84Ellipsoid


call geodetic2ecef(lat,lon,alt,x1,y1,z1)
call assert_isclose(x1,x0)
call assert_isclose(y1,y0)
call assert_isclose(z1,z0)

call aer2enu(az, el, rng, e1, n1, u1)
call assert_isclose(e1, er)
call assert_isclose(n1, nr)
call assert_isclose(u1, ur)

call aer2ecef(az,el,rng,lat,lon,alt,x2,y2,z2)
call assert_isclose(x2, xl)
call assert_isclose(y2, yl)
call assert_isclose(z2, zl)

call ecef2geodetic(x1,y1,z1,lat2,lon2,alt2)
call assert_isclose(lat2,lat)
call assert_isclose(lon2,lon)
call assert_isclose(alt2,alt)

call enu2aer(e1,n1,u1, az2, el2, rng2)
call assert_isclose(az2,az)
call assert_isclose(el2,el)
call assert_isclose(rng2,rng)

call ecef2aer(x2,y2,z2, lat,lon,alt, az3, el3, rng3)
call assert_isclose(az3, az)
call assert_isclose(el3, el)
call assert_isclose(rng3,rng)

call aer2geodetic(az,el,rng,lat,lon,alt, lat3,lon3,alt3)
call assert_isclose(lat3, lat1)

call geodetic2enu(lat3, lon3, alt3, lat, lon, alt, e2, n2, u2)
call assert_isclose(e2,e1)
call assert_isclose(n2,n1)
call assert_isclose(u2,u1)

call geodetic2aer(lat3,lon3,alt3,lat,lon,alt, az4, el4, rng4)
call assert_isclose(az4,az)
call assert_isclose(el4, el)
call assert_isclose(rng4, rng)

call enu2ecef(e1,n1,u1,lat,lon,alt, x3, y3, z3)
call assert_isclose(x3, x2)
call assert_isclose(y3, y2)
call assert_isclose(z3, z2)

call enu2geodetic(e2,n2,u2,lat,lon,alt,lat4, lon4, alt4)
call assert_isclose(lat4,lat3)
call assert_isclose(lon4,lon3)
call assert_isclose(alt4,alt3)

call ecef2enu(x3,y3,z3,lat,lon,alt, e3,n3,u3)
call assert_isclose(e3,e1)
call assert_isclose(n3,n1)
call assert_isclose(u3,u1)




end program
