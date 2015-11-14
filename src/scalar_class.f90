!=====================================================================!
! Module for switching between real and complex arithmentic modes
!=====================================================================!
module scalar_class

  ! module references
  use global_constants, only: dp

  ! module settings
  implicit none
  private 
  public :: scalar

  ! module definitions
  type scalar

     ! type settings
     private

     ! attributes
     real(dp)  :: x ! real part
     real(dp)  :: y ! imaginary part 

     ! procedures
   contains

     ! included below
     procedure :: get_real ! function
     procedure :: get_cplx ! function
     procedure :: set_real ! subroutine
     procedure :: set_cplx ! subroutine

  end type scalar

  ! interfaces
  interface scalar
     module procedure constructor
  end interface scalar

contains

  !-------------------------------------------------------------------!
  ! Constructor for scalar_class      
  !-------------------------------------------------------------------!
  function constructor(x, y) result(this)

    ! arguments
    real(dp)     :: x, y
    type(scalar) :: this

    call this%set_real(x)
    call this%set_cplx(y)

  end function constructor

  !-------------------------------------------------------------------!
  ! Getter for the real part
  !-------------------------------------------------------------------!
  function get_real(this)

    ! arguments
    real(dp)      :: get_real
    class(scalar) :: this

    ! return value is set
    get_real = this%x

  end function get_real

  !-------------------------------------------------------------------!
  ! getter for the complex part
  !-------------------------------------------------------------------!
  function get_cplx(this)

    ! arguments
    real(dp)      :: get_cplx
    class(scalar) :: this

    ! return value is set
    get_cplx = this%y

  end function get_cplx

  !-------------------------------------------------------------------!
  ! Setter for the real part
  !-------------------------------------------------------------------!
  subroutine set_real(this, x)

    ! arguments
    class(scalar) :: this
    real(dp)      :: x

    ! the value is set in the object
    this%x = x

  end subroutine set_real

  !-------------------------------------------------------------------!
  ! Setter for the real part
  !-------------------------------------------------------------------!
  subroutine set_cplx(this, y)

    ! arguments
    class(scalar) :: this
    real(dp)      :: y

    ! the value is set in the object
    this%y = y

  end subroutine set_cplx

end module scalar_class
