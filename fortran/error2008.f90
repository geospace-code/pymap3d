module error

implicit none

contains

pure subroutine errorstop
  error stop ! even Intel Fortran 2019 cannot handle string with error stop.
end subroutine errorstop

end module error

