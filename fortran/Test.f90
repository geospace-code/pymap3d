program Test

use iso_fortran_env, only: wp=>real64
use maptran

implicit none

real(wp), parameter :: lat = 42, lon= -82, alt = 200, &
                       x0 = 660.675e3, y0 = -4700.949e3, z0 = 4245.738e3

real(wp) :: lat2, lon2, alt2,x1,y1,z1

type(wgs84Ellipsoid) :: ell

ell = wgs84Ellipsoid()


call ecef2geodetic(x0,y0,z0,lat2,lon2,alt2)
call assert_isclose(lat2,lat)

call geodetic2ecef(lat,lon,alt,x1,y1,z1,ell,.true.)
call assert_isclose(x1,x0)
call assert_isclose(y1,y0)
call assert_isclose(z1,z0)


end program
