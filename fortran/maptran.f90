module maptran
  use iso_fortran_env, only: wp=>real64
  implicit none
  private

  type,public :: referenceEllipsoid
     real(wp) :: SemimajorAxis 
!     real :: Flattening = 1 / 298.2572235630  ! flattening
     real(wp) :: SemiminorAxis 
  end type

  real(wp), parameter :: pi = 4._wp * atan(1.0)
  
  type(referenceEllipsoid), parameter, public :: wgs84Ellipsoid = &
        referenceEllipsoid(SemimajorAxis=6378137., &
                           SemiminorAxis=6378137. * (1 - 1 / 298.2572235630))

  public :: ecef2geodetic, geodetic2ecef, aer2enu, enu2aer, aer2ecef, ecef2aer, &
            enu2ecef, ecef2enu, aer2geodetic, geodetic2enu,assert_isclose, &
            geodetic2aer,enu2geodetic

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
  type(referenceEllipsoid), intent(in), optional :: spheroid
  logical, optional,value :: deg
  real(wp), intent(out) :: lat, lon, alt

  real(wp) :: ea, eb, rad, rho, c, vnew, v
  integer :: i
  type(referenceEllipsoid) :: ell

  if (present(spheroid)) then
     ell = spheroid
  else
     ell = wgs84Ellipsoid
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
  type(referenceEllipsoid), intent(in), optional :: spheroid
  logical, optional, value :: deg

  real(wp) :: N
  type(referenceEllipsoid) :: ell

  if (present(spheroid)) then
     ell = spheroid
  else
     ell = wgs84Ellipsoid
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


subroutine aer2geodetic(az, el, slantRange, lat0, lon0, alt0, lat1, lon1, alt1, spheroid, deg)
!aer2geodetic  convert azimuth, elevation, range of target from observer to geodetic coordiantes
!
! Inputs
! ------
! az, el, slantrange: look angles and distance to point under test (degrees, degrees, meters)
! az: azimuth clockwise from local north
! el: elevation angle above local horizon
! lat0, lon0, alt0: ellipsoid geodetic coordinates of observer/reference (degrees, degrees, meters)
! spheroid: referenceEllipsoid parameter struct
! deg: .true.: degrees
!
! Outputs
! -------
! lat1,lon1,alt1: geodetic coordinates of test points (degrees,degrees,meters)
  
  real(wp), intent(in) :: az, el, slantRange, lat0, lon0, alt0 
  type(referenceEllipsoid), intent(in), optional :: spheroid
  logical, intent(in), optional :: deg
  real(wp), intent(out) :: lat1,lon1,alt1
  
  real(wp) :: x,y,z

  call aer2ecef(az, el, slantRange, lat0, lon0, alt0, x, y, z, spheroid, deg)
 
  call ecef2geodetic(x, y, z, lat1, lon1, alt1, spheroid, deg)
end subroutine aer2geodetic


subroutine geodetic2aer(lat, lon, alt, lat0, lon0, alt0, az, el, slantRange, spheroid, deg)
!geodetic2aer   from an observer's perspective, convert target coordinates to azimuth, elevation, slant range.
!
! Inputs
! ------
! lat,lon, alt:  ellipsoid geodetic coordinates of point under test (degrees, degrees, meters)
! lat0, lon0, alt0: ellipsoid geodetic coordinates of observer/reference (degrees, degrees, meters)
! spheroid: referenceEllipsoid parameter struct
! angleUnit: string for angular units. Default 'd': degrees, otherwise Radians
!
! Outputs
! -------
! az, el, slantrange: look angles and distance to point under test (degrees, degrees, meters)
! az: azimuth clockwise from local north
! el: elevation angle above local horizon

  real(wp), intent(in) :: lat,lon,alt, lat0, lon0, alt0 
  type(referenceEllipsoid), intent(in), optional :: spheroid
  logical, intent(in), optional :: deg
  real(wp), intent(out) :: az, el, slantRange

  real(wp) :: east,north,up


  call geodetic2enu(lat, lon, alt, lat0, lon0, alt0, east,north,up, spheroid, deg)
  call enu2aer(east, north, up, az, el, slantRange, deg)
  
end subroutine geodetic2aer



subroutine geodetic2enu(lat, lon, alt, lat0, lon0, alt0, east, north, up, spheroid, deg)
! geodetic2enu    convert from geodetic to ENU coordinates
!
! Inputs
! ------
! lat,lon, alt:  ellipsoid geodetic coordinates of point under test (degrees, degrees, meters)
! lat0, lon0, alt0: ellipsoid geodetic coordinates of observer/reference (degrees, degrees, meters)
! spheroid: referenceEllipsoid parameter struct
! angleUnit: string for angular units. Default 'd': degrees
!
! outputs
! -------
! e,n,u:  East, North, Up coordinates of test points (meters)

  real(wp), intent(in) :: lat, lon, alt, lat0, lon0, alt0 
  type(referenceEllipsoid), intent(in), optional :: spheroid
  logical, intent(in), optional :: deg
  real(wp), intent(out) :: east, north, up
  
  real(wp) x1,y1,z1,x2,y2,z2, dx,dy,dz


  call geodetic2ecef(lat,lon,alt,x1,y1,z1,spheroid,deg)
  call geodetic2ecef(lat0,lon0,alt0,x2,y2,z2,spheroid,deg)
  
  dx = x1-x2;
  dy = y1-y2;
  dz = z1-z2;
  
  call ecef2enuv(dx, dy, dz, lat0, lon0, east, north, up, deg)
  

end subroutine geodetic2enu


subroutine enu2geodetic(east, north, up, lat0, lon0, alt0, lat, lon, alt, spheroid, deg)
! enu2geodetic   convert from ENU to geodetic coordinates
!
! Inputs
! ------
!  East, North, Up: coordinates of point(s) (meters)
! lat0, lon0, alt0: ellipsoid geodetic coordinates of observer/reference (degrees, degrees, meters)
! spheroid: referenceEllipsoid parameter struct
! deg: .true. degrees
!
! outputs
! -------
! lat,lon,alt: geodetic coordinates of test points (degrees,degrees,meters)

  real(wp), intent(in) :: east, north, up, lat0, lon0, alt0 
  type(referenceEllipsoid), intent(in), optional :: spheroid
  logical, intent(in), optional :: deg
  real(wp), intent(out) :: lat, lon, alt
  
  real(wp) :: x,y,z

  call enu2ecef(east, north, up, lat0, lon0, alt0, x, y, z, spheroid, deg)
  call ecef2geodetic(x, y, z, lat, lon, alt, spheroid,deg)

end subroutine enu2geodetic



subroutine aer2ecef(az, el, slantRange, lat0, lon0, alt0, x,y,z, spheroid, deg)
! aer2ecef  convert azimuth, elevation, range to target from observer to ECEF coordinates
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
  type(referenceEllipsoid), intent(in), optional :: spheroid
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


subroutine ecef2aer(x, y, z, lat0, lon0, alt0, az, el, slantRange, spheroid, deg)
! ecef2aer  convert ECEF of target to azimuth, elevation, slant range from observer
!
! Inputs
! ------
! x,y,z: Earth Centered Earth Fixed (ECEF) coordinates of test point (meters)
! lat0, lon0, alt0: ellipsoid geodetic coordinates of observer/reference (degrees, degrees, meters)
! spheroid: referenceEllipsoid parameter struct
! deg: .true.: degrees
!
! Outputs
! -------
! az, el, slantrange: look angles and distance to point under test (degrees, degrees, meters)
! az: azimuth clockwise from local north
! el: elevation angle above local horizon

  real(wp), intent(in) :: x,y,z, lat0, lon0, alt0
  type(referenceEllipsoid), intent(in), optional :: spheroid
  logical, intent(in), optional :: deg
  real(wp), intent(out) :: az,el, slantRange
  
  real(wp) :: east, north, up

  call ecef2enu(x, y, z, lat0, lon0, alt0, east, north, up, spheroid, deg)
  call enu2aer(east, north, up, az, el, slantRange, deg)
  
end subroutine ecef2aer


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

  if (deg) az = radians(az); el = radians(el)

! Calculation of AER2ENU
   up = slantRange * sin(el)
   r = slantRange * cos(el)
   east = r * sin(az)
   north = r * cos(az)

end subroutine aer2enu


subroutine enu2aer(east, north, up, az, elev, slantRange, deg)
! enu2aer   convert ENU to azimuth, elevation, slant range
!
! Inputs
! ------
! e,n,u:  East, North, Up coordinates of test points (meters)
! deg: .true. degrees
!
! outputs
! -------
! az, el, slantrange: look angles and distance to point under test (degrees, degrees, meters)
! az: azimuth clockwise from local north
! el: elevation angle above local horizon

  real(wp),intent(in) :: east,north,up
  logical, optional, value :: deg
  real(wp), intent(out) :: az, elev, slantRange
  
  real(wp) :: r
  
  if (.not.present(deg)) deg = .true.

  r = hypot(east, north)
  slantRange = hypot(r, up)
  ! radians
  elev = atan2(up, r)
  az = mod(atan2(east, north), 2 * atan2(0.,-1.))

  if (deg) elev = degrees(elev); az = degrees(az)
  
end subroutine enu2aer


subroutine enu2ecef(e, n, u, lat0, lon0, alt0, x, y, z, spheroid, deg)
! enu2ecef  convert from ENU to ECEF coordiantes
!
! Inputs
! ------
! e,n,u:  East, North, Up coordinates of test points (meters)
! lat0, lon0, alt0: ellipsoid geodetic coordinates of observer/reference (degrees, degrees, meters)
! spheroid: referenceEllipsoid parameter struct
! angleUnit: string for angular units. Default 'd': degrees
!
! outputs
! -------
! x,y,z: Earth Centered Earth Fixed (ECEF) coordinates of test point (meters)
  real(wp), intent(in) :: e,n,u,lat0,lon0,alt0
  type(referenceEllipsoid), intent(in), optional :: spheroid
  logical, intent(in), optional :: deg
  real(wp), intent(out) :: x,y,z   
  
  real(wp) :: x0,y0,z0,dx,dy,dz           
            

  call geodetic2ecef(lat0, lon0, alt0, x0, y0, z0, spheroid, deg)
  call enu2uvw(e, n, u, lat0, lon0, dx, dy, dz, deg)
  
   x = x0 + dx
   y = y0 + dy
   z = z0 + dz
end subroutine enu2ecef



subroutine ecef2enu(x, y, z, lat0, lon0, alt0, east, north, up, spheroid, deg)
! ecef2enu  convert ECEF to ENU
!
! Inputs
! ------
! x,y,z: Earth Centered Earth Fixed (ECEF) coordinates of test point (meters)
! lat0, lon0, alt0: ellipsoid geodetic coordinates of observer/reference (degrees, degrees, meters)
! spheroid: referenceEllipsoid parameter struct
! angleUnit: string for angular units. Default 'd': degrees
!
! outputs
! -------
! e,n,u:  East, North, Up coordinates of test points (meters)
  
  real(wp), intent(in) :: x,y,z,lat0,lon0,alt0
  type(referenceEllipsoid), intent(in), optional :: spheroid
  logical, intent(in), optional :: deg
  real(wp), intent(out) :: east,north,up
  
  real(wp) :: x0,y0,z0

  call geodetic2ecef(lat0, lon0, alt0, x0,y0,z0,spheroid,deg)
  call ecef2enuv(x - x0, y - y0, z - z0, lat0, lon0, east, north, up, deg)
end subroutine ecef2enu


subroutine ecef2enuv(u, v, w, lat0, lon0, east, north, up, deg)
! ecef2enuv convert *vector projection* UVW to ENU
!
! Inputs
! ------
! u,v,w: meters
! lat0,lon0: geodetic latitude and longitude (degrees)
! deg: .true. degrees
!
! Outputs
! -------
! e,n,Up:  East, North, Up vector

  real(wp), intent(in) :: u,v,w
  real(wp), value :: lat0,lon0
  logical, optional, value :: deg
  real(wp), intent(out) :: east, north, up
  
  real(wp) :: t
  
  if (.not.present(deg)) deg = .true.
  
  if (deg) lat0 = radians(lat0); lon0 = radians(lon0)
  
  t     =  cos(lon0) * u + sin(lon0) * v
  east  = -sin(lon0) * u + cos(lon0) * v
  up    =  cos(lat0) * t + sin(lat0) * w
  north = -sin(lat0) * t + cos(lat0) * w
end subroutine ecef2enuv


subroutine enu2uvw(e,n,up,lat0,lon0,u,v,w,deg)
! enu2uvw   convert from ENU to UVW coordinates
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
    type(referenceEllipsoid), intent(in) :: E

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
