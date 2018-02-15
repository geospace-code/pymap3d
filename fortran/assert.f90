module assert

  use, intrinsic:: iso_fortran_env, only: sp=>real32, dp=>real64, qp=>real128, stderr=>error_unit
  use, intrinsic:: ieee_arithmetic
  implicit none
  private
  
  integer,parameter :: wp = dp
  
  public :: wp,isclose, assert_isclose
  
contains

elemental logical function isclose(actual, desired, rtol, atol, equal_nan)
! inputs
! ------
! actual: value "measured"
! desired: value "wanted"
! rtol: relative tolerance
! atol: absolute tolerance
! equal_nan: consider NaN to be equal?
!
!  rtol overrides atol when both are specified
!
! https://www.python.org/dev/peps/pep-0485/#proposed-implementation
! https://github.com/PythonCHB/close_pep/blob/master/is_close.py

  real(wp), intent(in) :: actual, desired
  real(wp), intent(in), optional :: rtol, atol
  logical, intent(in), optional :: equal_nan
  
  real(wp) :: r,a
  logical :: n

  r = merge(rtol, 1e-5_wp, present(rtol)) 
  a = merge(atol, 0._wp, present(atol)) 
  n = merge(equal_nan, .false., present(equal_nan))
  
  !print*,r,a,n,actual,desired
  
!--- sanity check
  if ((r < 0._wp).or.(a < 0._wp)) error stop 'tolerances must be non-negative'
!--- simplest case
  isclose = (actual == desired) 
  if (isclose) return
!--- equal nan
  isclose = n.and.(ieee_is_nan(actual).and.ieee_is_nan(desired))
  if (isclose) return
!--- Inf /= -Inf, unequal NaN
  if (.not.ieee_is_finite(actual) .or. .not.ieee_is_finite(desired)) return
!--- floating point closeness check
  isclose = abs(actual-desired) <= max(r * max(abs(actual), abs(desired)), a)

end function isclose


impure elemental subroutine assert_isclose(actual, desired, rtol, atol, equal_nan)
! inputs
! ------
! actual: value "measured"
! desired: value "wanted"
! rtol: relative tolerance
! atol: absolute tolerance
! equal_nan: consider NaN to be equal?
!
! rtol overrides atol when both are specified

  real(wp), intent(in) :: actual, desired
  real(wp), intent(in), optional :: rtol, atol
  logical, intent(in), optional :: equal_nan
 
  if (.not.isclose(actual,desired,rtol,atol,equal_nan)) then
    write(stderr,*) 'actual',actual,'desired',desired
    error stop
  endif

end subroutine assert_isclose

end module assert
