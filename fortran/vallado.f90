module vallado
! based on http://www.smad.com/vallado/

use assert, only : wp

implicit none
private

  real(wp), parameter :: pi = 4._wp * atan(1.0_wp)

  type,public :: time
     integer :: year, month, day, hour, minute
     real(wp) :: second
  end type

public:: toGST, toJulian, toLST, radec2azel, azel2radec
contains

elemental subroutine azel2radec(az,el,lat,lon,jd, ra,decl)
! convert azimuth, elevation to right ascension, declination
! 
! inputs
! ------
! az, el: azimuth, elevation (degrees)
! lat, lon: geodetic latitude, longitude (degrees)
! jd: julian date (decimal)
! 
! outputs
! -------
! ra, decl: right ascension, declination (degrees)

  real(wp), value :: az,el,lat,lon
  real(wp), intent(in) :: jd ! Julian Date
  real(wp), intent(out) :: ra, decl
  
  real(wp) :: lst, lha, sinv, cosv
  
  az = radians(az)
  el = radians(el)
  lat = radians(lat)
  lon = radians(lon)

  Decl = ASIN(SIN(El)*SIN(lat) + &
              COS(el)*COS(lat)*COS(Az) )

  Sinv = -(SIN(az)*COS(el)*COS(lat)) / &
          (COS(lat)*COS(Decl))
          
  Cosv = (SIN(el) - SIN(lat)*SIN(decl)) / &
         (COS(lat)*COS(Decl))
         
  LHA  = ATAN2(Sinv,Cosv)
  lst = toLST(Lon, JD)
  
  ra = modulo(degrees(LST - LHA), 360._wp)
  decl = degrees(decl)

end subroutine azel2radec


elemental SUBROUTINE radec2azel(ra,Decl,lat,lon,jd, Az,El)
! convert right ascension, declination to azimuth, elevation 
! 
! inputs
! ------
! ra, decl: right ascension, declination (degrees)
! lat, lon: geodetic latitude, longitude (degrees)
! jd: julian date (decimal)
! 
! outputs
! -------
! az, el: azimuth, elevation (degrees)

  REAL(wp), value :: ra,Decl,lat, lon
  real(wp), intent(in) :: jd
  real(wp), intent(out) :: Az,El
  REAL(wp) :: Sinv, Cosv, LHA
  
  lat = radians(lat)
  lon = radians(lon)
  ra = radians(ra)
  decl = radians(decl)

  LHA = toLST(Lon, JD) - ra

  El = ASIN( SIN(Decl)*SIN(lat) + &
              COS(Decl)*COS(lat)*COS(LHA) )

  Sinv = -SIN(LHA)*COS(Decl)*COS(lat) / &
        (COS(el)*COS(lat))
        
  Cosv = ( SIN(Decl)-SIN(el)*SIN(lat) ) / &
         (COS(el)*COS(lat))
         
  Az = modulo(degrees(ATAN2(Sinv,Cosv)), 360._wp)
  el = degrees(el)


END subroutine radec2azel


elemental real(wp) function toLST(Lon, JD) result(LST)
! Julian Date => local sidereal time
!
! inputs
! ------
! lon: geodetic longitude (radians)
! jd: Julian Date (decimal)

  REAL(wp), intent(in) :: Lon, JD
 
  LST = Lon + toGST(jd)

  LST = modulo(LST, 2*pi )

END function toLST
      

elemental real(wp) function toJulian(t) result(jd)
! Gregorian date, time => Julian Date
!
! inputs
! ------
! time: Gregorian user-defined type
!
! output
! -----
! JD: Julian Date

  type(time), intent(in) :: t
  real(wp) :: B, y, m
  
  y = t%year
  m = t%month

  IF ( M <= 2 ) THEN
    Y = y - 1
    M = m + 12
  ENDIF
    
  B = 2 - INT(Y*0.01_wp) + INT(INT(Y*0.01_wp)*0.25_wp)
  
  JD= INT( 365.25_wp*(Y + 4716) ) + &
      INT( 30.6001_wp*(M+1) ) + &
      t%day + B - 1524.5_wp + &
      ( (t%second/60.0_wp + t%minute ) / 60.0_wp + t%hour ) / 24.0_wp

END function toJulian


elemental real(wp) FUNCTION toGST(JD) result(GST)
! Julian Date => to Greenwich Sidereal Time
!
! inputs
! ------
! JD: Julian Date (decimal)
!
! output
! ------
! GST: Greenwich Sidereal Time (decimal)

  real(wp), intent(in) :: JD
  real(wp) :: TUT1


  TUT1= ( JD - 2451545._wp ) / 36525._wp
  gst = -6.2e-6_wp*TUT1**3 + &
              0.093104_wp*TUT1**2 + &
             (876600._wp*3600._wp + 8640184.812866_wp)*TUT1 + &
              67310.54841_wp
              
  gst = modulo(radians(gst) / 240._wp, 2*pi) ! 360/86400 = 1/240, to deg, to rad
  

end function toGST


elemental real(wp) function degrees(rad)
    real(wp), intent(in) :: rad

    degrees = 180._wp / pi * rad
end function degrees


elemental real(wp) function radians(deg)
    real(wp), intent(in) :: deg

    radians = pi / 180._wp * deg
end function radians


end module
