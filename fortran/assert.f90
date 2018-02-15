module assert

  use, intrinsic:: iso_fortran_env, only: wp=>real64, stderr=>error_unit
  implicit none
  private
  
  public :: isclose, assert_isclose
  
contains

elemental logical function isclose(actual, desired, rtol, atol)
! https://www.python.org/dev/peps/pep-0485/#proposed-implementation
! Value wasn't used for atol and rtol to keep this fundamental function "pure"
    real(wp), intent(in) :: actual, desired
    real(wp), optional, value:: rtol, atol
    real(wp) :: r,a

    if (.not.present(rtol)) then
        r = 1e-5_wp
    else
        r = rtol
    endif
    
    if (.not.present(atol)) then
        a = 0._wp
    else
        a = atol
    endif

    isclose = (abs(actual-desired) <= max(r * max(abs(actual), abs(desired)), a))
end function isclose


impure elemental subroutine assert_isclose(actual, desired, rtol, atol)
    real(wp), intent(in) :: actual, desired
    real(wp), intent(in), optional :: rtol, atol
    logical ok
    ok = isclose(actual,desired,rtol,atol)

    if (.not.ok) then
        write(stderr,*) 'actual',actual,'desired',desired
        error stop
    endif

end subroutine assert_isclose

end module assert
