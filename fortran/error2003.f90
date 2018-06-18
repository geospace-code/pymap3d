module error

implicit none

public :: errorstop

contains

  pure subroutine errorstop
  
    stop 1
  end subroutine errorstop

end module error

