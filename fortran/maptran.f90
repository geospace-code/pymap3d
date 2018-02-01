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

  public :: ecef2geodetic, geodetic2ecef, aer2enu, aer2ecef, assert_isclose

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
  logical, optional,value :: deg
  real(wp), intent(out) :: lat, lon, alt

  real(wp) :: ea, eb, rad, rho, c, vnew, v
  integer :: i
  type(wgs84Ellipsoid) :: ell

  if (present(spheroid)) then
     ell = spheroid
  else
     ell = wgs84Ellipsoid()
  endif

  if (.not.present(deg)) deg = .true.

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

    if (deg) lat = degrees(lat); lon = degrees(lon)

end subroutine ecef2geodetic


subroutine geodetic2ecef(lat,lon,alt,x,y,z,spheroid,deg)
  real(wp), value :: lat,lon
  real(wp), intent(in) :: alt
  real(wp), intent(out) :: x,y,z
  type(wgs84Ellipsoid), intent(in), optional :: spheroid
  logical, optional, value :: deg

  real(wp) :: N
  type(wgs84Ellipsoid) :: ell

  if (present(spheroid)) then
     ell = spheroid
  else
     ell = wgs84Ellipsoid()
  endif

  if (.not.present(deg)) deg=.true.

  if (deg) lat = radians(lat); lon = radians(lon)

! Radius of curvature of the prime vertical section
  N = radius_normal(lat, ell)
! Compute cartesian (geocentric) coordinates given  (curvilinear) geodetic coordinates.

  x = (N + alt) * cos(lat) * cos(lon)
  y = (N + alt) * cos(lat) * sin(lon)
  z = (N * (ell%SemiminorAxis / ell%SemimajorAxis)**2 + alt) * sin(lat)

end subroutine geodetic2ecef


subroutine aer2ecef(az, el, slantRange, lat0, lon0, alt0, x,y,z, spheroid, deg)
!er2ecef  convert azimuth, elevation, range to target from observer to ECEF coordinates
!
! Inputs
! ------
! az, el, slantrange: look angles and distance to point under test (degrees, degrees, meters)
! az: azimuth clockwise from local north
! el: elevation angle above local horizon
! lat0, lon0, alt0: ellipsoid geodetic coordinates of observer/reference (degrees, degrees, meters)
! spheroid: referenceEllipsoid parameter struct
! deg: .true. degrees
!
! outputs
! -------
! x,y,z: Earth Centered Earth Fixed (ECEF) coordinates of test point (meters)

  real(wp), intent(in) :: az,el, slantRange, lat0, lon0, alt0
  type(wgs84Ellipsoid), intent(in), optional :: spheroid
  logical, intent(in), optional :: deg
  real(wp), intent(out) :: x,y,z

  real(wp) :: x0,y0,z0, e,n,u,dx,dy,dz


! Origin of the local system in geocentric coordinates.
  call geodetic2ecef(lat0, lon0, alt0,x0, y0, z0, spheroid,deg)
! Convert Local Spherical AER to ENU
  call aer2enu(az, el, slantRange, e, n, u,deg)
! Rotating ENU to ECEF
  call enu2uvw(e, n, u, lat0, lon0, dx, dy, dz,deg)
! Origin + offset from origin equals position in ECEF
  x = x0 + dx
  y = y0 + dy
  z = z0 + dz

end subroutine aer2ecef


subroutine aer2enu(az, el, slantRange, east, north, up, deg)
!aer2enu  convert azimuth, elevation, range to ENU coordinates
!
! Inputs
! ------
! az, el, slantrange: look angles and distance to point under test (degrees, degrees, meters)
! az: azimuth clockwise from local north
! el: elevation angle above local horizon
! angleUnit: string for angular units. Default 'd': degrees
!
! Outputs
! -------
! e,n,u:  East, North, Up coordinates of test points (meters)

  real(wp), value :: az,el
  real(wp), intent(in) :: slantRange
  logical, optional, value :: deg
  real(wp),intent(out) :: east, north,up

  real(wp) :: r

  if(.not.present(deg)) deg=.true.

  if (deg) az = degrees(az); el = degrees(el)

! Calculation of AER2ENU
   up = slantRange * sin(el)
   r = slantRange * cos(el)
   east = r * sin(az)
   north = r * cos(az)

end subroutine aer2enu


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
  real(wp), intent(in) :: e,n,up
  real(wp), value :: lat0,lon0
  real(wp), intent(out) :: u,v,w
  logical, optional, value :: deg

  real(wp) :: t

  if(.not.present(deg)) deg=.true.

  if (deg) lat0 = radians(lat0); lon0 = radians(lon0)


    t = cos(lat0) * up - sin(lat0) * n
    w = sin(lat0) * up + cos(lat0) * n

    u = cos(lon0) * t - sin(lon0) * e
    v = sin(lon0) * t + cos(lon0) * e

end subroutine enu2uvw


elemental real(wp) function radius_normal(lat,E)
    real(wp), intent(in) :: lat
    type(wgs84Ellipsoid), intent(in) :: E

    radius_normal = E%SemimajorAxis**2 / sqrt( E%SemimajorAxis**2 * cos(lat)**2 + E%SemiminorAxis**2 * sin(lat)**2 )

end function radius_normal


logical function isclose(actual, desired, rtol, atol)
! https://www.python.org/dev/peps/pep-0485/#proposed-implementation
    real(wp), intent(in) :: actual, desired
    real(wp), optional,value :: rtol, atol

    if (.not.present(rtol)) rtol = 0.01
    if (.not.present(atol)) atol = 0.

    isclose = (abs(actual-desired) <= max(rtol * max(abs(actual), abs(desired)), atol))
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
