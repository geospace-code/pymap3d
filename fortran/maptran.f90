module maptran
  use iso_fortran_env, only: wp=>real64
  implicit none
  private

  type, public :: wgs84Ellipsoid
     real(wp) :: SemimajorAxis = 6378137.  ! [m]
!     real :: Flattening = 1 / 298.2572235630  ! flattening
     real(wp) :: SemiminorAxis = 6378137. * (1 - 1 / 298.2572235630)
  end type

  real(wp), parameter :: pi = 4._wp * atan(1.0)
  
  public :: ecef2geodetic, geodetic2ecef, assert_isclose

contains

subroutine ecef2geodetic(x, y, z, lat, lon, alt, spheroid, deg)

! convert ECEF (meters) to geodetic coordintes

! Algorithm is based on
! http://www.astro.uni.torun.pl/~kb/Papers/geod/Geod-BG.htm
! This algorithm provides a converging solution to the latitude equation
! in terms of the parametric or reduced latitude form (v)
! This algorithm provides a uniform solution over all latitudes as it does
! not involve division by cos(phi) or sin(phi)

  real(wp), intent(in) :: x,y,z
  type(wgs84Ellipsoid), intent(in), optional :: spheroid
  logical, intent(in), optional :: deg
  real(wp), intent(out) :: lat, lon, alt

  type(wgs84Ellipsoid) :: ell
  logical dg
  real(wp) :: ea, eb, rad, rho, c, vnew, v
  integer :: i

  if (present(spheroid)) then
     ell = spheroid
  else
     ell = wgs84Ellipsoid()
  endif

  if (present(deg)) then
      dg = deg
  else
      dg = .true.
  endif

    ea = ell%SemimajorAxis
    eb = ell%SemiminorAxis
    rad = hypot(x, y)
! Constant required for Latitude equation
    rho = atan2(eb * z, ea * rad)
! Constant required for latitude equation
    c = (ea**2 - eb**2) / hypot(ea * rad, eb * z)
! Starter for the Newtons Iteration Method
    vnew = atan2(ea * z, eb * rad)
! Initializing the parametric latitude
    v = 0.
    do i = 1,5
      v = vnew
     ! Newtons Method for computing iterations
      vnew = v - ((2 * sin(v - rho) - c * sin(2 * v)) / &
                (2 * (cos(v - rho) - c * cos(2 * v))))

      if (isclose(v,vnew)) exit
    enddo

! Computing latitude from the root of the latitude equation
    lat = atan2(ea * tan(vnew), eb)
! by inspection
    lon = atan2(y, x)

    alt = ((rad - ea * cos(vnew)) * cos(lat)) + &
           ((z - eb * sin(vnew)) * sin(lat))

    if (dg) lat = degrees(lat); lon = degrees(lon)

end subroutine ecef2geodetic


subroutine geodetic2ecef(lat,lon,alt,x,y,z,spheroid,deg)
  real(wp), intent(in) :: lat,lon,alt
  real(wp), intent(out) :: x,y,z
  type(wgs84Ellipsoid) :: spheroid
  logical, intent(in) :: deg
  
  real(wp) :: N, lat1,lon1
  
  if (deg) then
    lat1 = radians(lat); lon1 = radians(lon)
  else
    lat1 = lat; lon1 = lon
  endif
! Radius of curvature of the prime vertical section
  N = radius_normal(lat1, spheroid)
! Compute cartesian (geocentric) coordinates given  (curvilinear) geodetic coordinates.
  
  x = (N + alt) * cos(lat1) * cos(lon1)
  y = (N + alt) * cos(lat1) * sin(lon1)
  z = (N * (spheroid%SemiminorAxis / spheroid%SemimajorAxis)**2 + alt) * sin(lat1)

end subroutine geodetic2ecef


subroutine enu2uvw(e,n,up,lat0,lon0,u,v,w,deg)
!enu2uvw   convert from ENU to UVW coordinates
!
! Inputs
! ------
! e,n,up:  East, North, Up coordinates of point(s) (meters)
! lat0,lon0: geodetic coordinates of observer/reference point (degrees)
! deg: ,true. degrees
!
! outputs
! -------
! u,v,w:   coordinates of test point(s) (meters)
  real(wp), intent(in) :: e,n,up,lat0,lon0
  real(wp), intent(out) :: u,v,w
  logical, intent(in) :: deg
  
  real(wp) :: t,lat,lon
  
  if (deg) then
    lat = radians(lat0); lon = radians(lon0)
  else
    lat = lat0; lon= lon0
  endif

    
    t = cos(lat) * up - sin(lat) * n
    w = sin(lat) * up + cos(lat) * n

    u = cos(lon) * t - sin(lon) * e
    v = sin(lon) * t + cos(lon) * e
    
end subroutine enu2uvw


elemental real(wp) function radius_normal(lat,E)
    real(wp), intent(in) :: lat
    type(wgs84Ellipsoid), intent(in) :: E

    radius_normal = E%SemimajorAxis**2 / sqrt( E%SemimajorAxis**2 * cos(lat)**2 + E%SemiminorAxis**2 * sin(lat)**2 )

end function radius_normal 


elemental logical function isclose(actual, desired, rtol, atol)
! https://www.python.org/dev/peps/pep-0485/#proposed-implementation
    real(wp), intent(in) :: actual, desired
    real(wp), intent(in), optional :: rtol, atol
    real(wp) :: r, a

    if (present(rtol)) then
        r = rtol
    else
        r = 1.0e-6
    endif

    if (present(atol)) then
        a = atol
    else
        a = 0.
    endif

    isclose = (abs(actual-desired) <= max(r * max(abs(actual), abs(desired)), a))
end function isclose


subroutine assert_isclose(actual, desired, rtol, atol)
    real(wp), intent(in) :: actual, desired
    real(wp), intent(in), optional :: rtol, atol
    logical ok
    ok = isclose(actual,desired,rtol,atol)
    
    if (.not.ok) then
        print*,'actual',actual,'desired',desired
        error stop
    endif
    
end subroutine assert_isclose


elemental real(wp) function degrees(rad)
    real(wp), intent(in) :: rad

    degrees = 180._wp / pi * rad
end function degrees


elemental real(wp) function radians(deg)
    real(wp), intent(in) :: deg

    radians = pi / 180._wp * deg
end function radians

endmodule
