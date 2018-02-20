module vallado

use, intrinsic:: iso_fortran_env, only: dp=>real64, sp=>real32
implicit none
private
integer, parameter :: wp=dp
real(wp), parameter :: pi = 4._wp * atan(1.0_wp)

public:: gstime, juliandayall, lstime
contains


elemental SUBROUTINE RADEC_AZEL  ( RtAsc,Decl,LST,Latgd, Direction, Az,El )

  REAL(wp), value :: RtAsc,Decl,LST,Latgd
  real(wp), intent(out) :: Az,El
  CHARACTER(4), intent(in), optional :: Direction
  REAL(wp) :: Sinv, Cosv, LHA

        IF (present(direction) .and. Direction .eq. 'FROM' ) THEN
            Decl = ASIN( SIN(El)*SIN(Latgd) + &
     &                 COS(el)*COS(Latgd)*COS(Az) )

            Sinv = -(SIN(az)*COS(el)*COS(Latgd)) / &
     &              (COS(Latgd)*COS(Decl))
            Cosv = (SIN(el) - SIN(Latgd)*SIN(decl)) / &
     &              (COS(Latgd)*COS(Decl))
            LHA  = ATAN2( Sinv,Cosv ) 
            RtAsc= LST - LHA 
          ELSE
            LHA = LST - RtAsc

            El  = ASIN( SIN(Decl)*SIN(Latgd) + &
     &            COS(Decl)*COS(Latgd)*COS(LHA) )

            Sinv= -SIN(LHA)*COS(Decl)*COS(Latgd) / &
     &                (COS(el)*COS(Latgd))
            Cosv= ( SIN(Decl)-SIN(el)*SIN(Latgd) ) / &
     &             (COS(el)*COS(Latgd))
            Az  = ATAN2( Sinv,Cosv ) 
          ENDIF 


END subroutine radec_azel


elemental real(wp) function LSTIME( Lon, JD)

  REAL(wp), intent(in) :: Lon, JD
  real(wp) :: GST


  GST = GSTIME(jd)
  LSTime = Lon + GST

  LSTime = MOD( LSTime, 2*pi )
  IF ( LSTime < 0._wp ) LSTime = LSTime + 2*pi

END function lstime
      

elemental real(wp) function JULIANDAYALL(Year,Mon,Day,Hr,Mint,Sec, WhichType) result(jd)

  INTEGER,intent(in) :: Year, Mon, Day, Hr, Mint
  real(wp), intent(in) :: sec
  CHARACTER, intent(in), optional :: WhichType
  real(wp) :: B, y, m
  
  y = year
  m = mon

  IF ( M <= 2 ) THEN
      Y = y - 1
      M = m + 12
    ENDIF
  IF (present(whichtype).and. WhichType == 'J') THEN
      ! -------- Use for Julian calender, every 4 years ---------
      B = 0._wp
  ELSE
      ! --------------------- Use for Gregorian -----------------
      B = 2 - INT(Y*0.01_wp) + INT(INT(Y*0.01_wp)*0.25_wp)
  ENDIF
  
  JD= INT( 365.25_wp*(Y + 4716) ) + &
      INT( 30.6001_wp*(M+1) ) + &
      Day + B - 1524.5_wp + &
      ( (Sec/60.0_wp + Mint ) / 60.0_wp + Hr ) / 24.0_wp

END function juliandayall


pure elemental real(wp) FUNCTION GSTIME(JD)
! julian2sidereal(100000) = 2.9310980581630943


  real(wp), intent(in) :: JD
  real(wp) :: TUT1


  TUT1= ( JD - 2451545._wp ) / 36525._wp
  gstime = -6.2e-6_wp*TUT1**3 + &
              0.093104_wp*TUT1**2 + &
             (876600._wp*3600._wp + 8640184.812866_wp)*TUT1 + &
              67310.54841_wp
              
  gstime = MOD(radians(gstime) / 240._wp, 2*pi) ! 360/86400 = 1/240, to deg, to rad
  
  IF (gstime < 0._wp) gstime = gstime + 2*pi


end function gstime


elemental real(wp) function degrees(rad)
    real(wp), intent(in) :: rad

    degrees = 180._wp / pi * rad
end function degrees


elemental real(wp) function radians(deg)
    real(wp), intent(in) :: deg

    radians = pi / 180._wp * deg
end function radians


end module
