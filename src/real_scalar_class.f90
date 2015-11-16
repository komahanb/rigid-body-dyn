!=====================================================================!
! Module for switching between real and complex arithmentic modes
!=====================================================================!
module real_scalar_class

  ! module references
  use global_constants, only : dp
  use scalar_class, only : scalar

  ! module settings
  implicit none
  private 
  public :: real_scalar

  ! module definition
  type, extends(scalar) :: real_scalar

     ! nothing changed from the parent
  end type real_scalar
  
  !interface
  interface real_scalar
     procedure constructor
  end interface real_scalar
  
  ! module procedure
contains
  
  !-------------------------------------------------------------------!
  ! constructor for real_scalar_class      
  !-------------------------------------------------------------------!
  function constructor(x) result(this)

    ! arguments
    real(dp)          :: x
    type(real_scalar) :: this
    
    call this%set_real(x)
    
  end function constructor

end module real_scalar_class
