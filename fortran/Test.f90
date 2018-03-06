program Test

use maptran
use vallado
use assert, only:wp, assert_isclose

implicit none

!type(Ellipsoid), parameter :: spheroid = wgs84Ellipsoid

real(wp), parameter :: lat = 42, lon= -82, alt = 200, &
                       az = 33, el=70, rng= 1e3_wp, &
                       x0 = 660.675e3_wp, y0 = -4700.949e3_wp, z0 = 4245.738e3_wp, &
                       xl = 660.930e3_wp, yl = -4701.424e3_wp, zl = 4246.579e3_wp, & ! aer2ecef
                       er = 186.277521_wp, nr = 286.84222_wp, ur = 939.69262_wp, & ! aer2enu
                       lat1 = 42.0026_wp, lon1 = -81.9978_wp, alt1 = 1.1397e3_wp, & ! aer2geodetic
                       azi = 180.1_wp,  eli = 80._wp, &
                       ha = 45.482789587392013_wp
               
integer,parameter :: N = 3               
real(wp), dimension(N), parameter :: alat = [42,52,62], &
                                     deg0 = [15,30,45], &
                                     aaz = [33,43,53]
                                     
type(time), parameter :: t0 = time(2014,4,6,8,0,0) ! UTC 
real(wp), parameter ::  jd0 = 2456753.833333_wp



real(wp) :: lat2, lon2, alt2,lat3,lon3,alt3,lat4,lon4,alt4,&
            x1,y1,z1,x2,y2,z2,x3,y3,z3,&
            az2,el2,rng2,az3,el3,rng3,az4,el4,rng4,azrd,elrd,rae,dae,jd, &
            e1,n1,u1,e2,n2,u2,e3,n3,u3
  

real(wp), dimension(N) :: ax1, ay1, aaaz1, ax2, ay2, aaaz2, ax3,ay3,aaaz3, &
                          ae1, an1, au1, ae2,an2,au2, &
                          alat2,alon2,aalt2, alat3,alon3,aalt3, alat4,alon4,aalt4, &
                          aaz2, ael2, arng2, aaz3,ael3,arng3, aaz4,ael4,arng4
                         


!print*,'Default WGS84 Ellipsoid:',spheroid

! --------- scalar degrees

call geodetic2ecef(lat,lon,alt, x1,y1,z1)
call assert_isclose([x1,y1,z1],[x0,y0,z0])

call aer2enu(az,el,rng, e1,n1,u1)
call assert_isclose([e1,n1,u1], [er,nr,ur])

call aer2ecef(az,el,rng,lat,lon,alt,x2,y2,z2)
call assert_isclose([x2,y2,z2],[xl,yl,zl])

call ecef2geodetic(x1,y1,z1,lat2,lon2,alt2)
call assert_isclose([lat2,lon2,alt2],[lat,lon,alt])

call enu2aer(e1,n1,u1, az2, el2, rng2)
call assert_isclose([az2,el2,rng2],[az,el,rng])

call ecef2aer(x2,y2,z2, lat,lon,alt, az3, el3, rng3)
call assert_isclose([az3,el3,rng3],[az,el,rng])

call aer2geodetic(az,el,rng,lat,lon,alt, lat3,lon3,alt3)
call assert_isclose([lat3,lon3,alt3],[lat1,lon1,alt1])

call geodetic2enu(lat3, lon3, alt3, lat, lon, alt, e2, n2, u2)
call assert_isclose([e2,n2,u2],[e1,n1,u1])

call geodetic2aer(lat3,lon3,alt3,lat,lon,alt, az4, el4, rng4)
call assert_isclose([az4,el4,rng4],[az,el,rng])

call enu2ecef(e1,n1,u1,lat,lon,alt, x3, y3, z3)
call assert_isclose([x3,y3,z3],[x2,y2,z2])

call enu2geodetic(e2,n2,u2,lat,lon,alt,lat4, lon4, alt4)
call assert_isclose([lat4,lon4,alt4],[lat3,lon3,alt3])

call ecef2enu(x3,y3,z3,lat,lon,alt, e3,n3,u3)
call assert_isclose([e3,n3,u3],[e1,n1,u1])

call assert_isclose(degrees(radians(deg0)), deg0)

! ------ scalar radians

call geodetic2ecef(radians(lat),radians(lon),alt, x1,y1,z1, deg=.false.)
call assert_isclose([x1,y1,z1],[x0,y0,z0])

call aer2enu(radians(az),radians(el),rng, e1,n1,u1, deg=.false.)
call assert_isclose([e1,n1,u1], [er,nr,ur])

call aer2ecef(radians(az),radians(el),rng, radians(lat),radians(lon),alt, x2,y2,z2, deg=.false.)
call assert_isclose([x2,y2,z2],[xl,yl,zl])

call ecef2geodetic(x1,y1,z1,lat2,lon2,alt2, deg=.false.)
call assert_isclose([degrees(lat2),degrees(lon2),alt2],[lat,lon,alt])

call enu2aer(e1,n1,u1, az2, el2, rng2, deg=.false.)
call assert_isclose([degrees(az2),degrees(el2),rng2],[az,el,rng], err_msg='enu2aer: rad')

call ecef2aer(x2,y2,z2, radians(lat),radians(lon),alt, az3,el3,rng3, deg=.false.)
call assert_isclose([degrees(az3),degrees(el3),rng3],[az,el,rng])

call aer2geodetic(radians(az),radians(el),rng, radians(lat),radians(lon),alt, lat3,lon3,alt3, deg=.false.)
call assert_isclose([degrees(lat3),degrees(lon3),alt3],[lat1,lon1,alt1], err_msg='aer2geodetic: rad')

call geodetic2enu(lat3,lon3,alt3, radians(lat),radians(lon),alt, e2,n2,u2, deg=.false.)
call assert_isclose([e2,n2,u2],[e1,n1,u1], err_msg='geodetic2enu: rad')

call geodetic2aer(lat3,lon3,alt3, radians(lat),radians(lon),alt, az4,el4,rng4, deg=.false.)
call assert_isclose([degrees(az4),degrees(el4),rng4],[az,el,rng], err_msg='geodetic2aer: rad')

call enu2ecef(e1,n1,u1, radians(lat),radians(lon),alt, x3,y3,z3, deg=.false.)
call assert_isclose([x3,y3,z3],[x2,y2,z2], err_msg='enu2ecef: rad')

call enu2geodetic(e2,n2,u2, radians(lat),radians(lon),alt, lat4,lon4,alt4, deg=.false.)
call assert_isclose([lat4,lon4,alt4],[lat3,lon3,alt3], err_msg='enu2geodetic: rad')

call ecef2enu(x3,y3,z3, radians(lat),radians(lon),alt, e3,n3,u3, deg=.false.)
call assert_isclose([e3,n3,u3],[e1,n1,u1], err_msg='ecef2enu: rad')

! --- array



call geodetic2ecef(alat,lon,alt,ax1,ay1,aaaz1)
call aer2enu(aaz, el, rng, ae1, an1, au1)
call aer2ecef(aaz, el, rng, lat,lon,alt, ax2, ay2, aaaz2)
call ecef2geodetic(ax1,ay1,aaaz1,alat2,alon2,aalt2)
call enu2aer(ae1,an1,au1, aaz2, ael2, arng2)
call ecef2aer(ax2,ay2,az2, lat,lon,alt, aaz3, ael3, arng3)
call aer2geodetic(aaz,el,rng,lat,lon,alt, alat3,alon3,aalt3)
call geodetic2enu(alat3, alon3, aalt3, lat, lon, alt, ae2, an2, au2)
call geodetic2aer(alat3,alon3, aalt3,lat,lon,alt, aaz4, ael4, arng4)
call enu2ecef(ae1,an1,au1,lat,lon,alt, ax3, ay3, aaaz3)
call enu2geodetic(ae2,an2,au2,lat,lon,alt,alat4, alon4, aalt4)

!-------- Vallado
jd = toJulian(t0)
! http://aa.usno.navy.mil/jdconverter?ID=AA&year=2014&month=4&day=6&era=1&hr=8&min=0&sec=0.0
call assert_isclose(jd, jd0, err_msg='toJulian')

call assert_isclose(toGST(jd), 5.4896448816_wp, err_msg='toGST')
call assert_isclose(toLST(radians(-148._wp), jd), 2.90658_wp, err_msg='toLST')

call azel2radec(azi,eli,lat,lon, jd, rae, dae)
call radec2azel(rae,dae,lat,lon,jd,azrd,elrd)
call assert_isclose([azi,eli],[azrd,elrd])
!------- Meeus
call assert_isclose(anglesep(35._wp,23._wp, 84._wp,20._wp), ha)


print *,'Maptran OK'
end program
