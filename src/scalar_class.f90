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
     real(dp)  :: x = 0.0_dp ! real part
     real(dp)  :: y = 0.0_dp ! imaginary part 

     ! type procedures
   contains

     private

     ! included below
     procedure, public :: get_real ! function
     procedure, public :: get_cplx ! function
     procedure, public :: set_real ! subroutine
     procedure, public :: set_cplx ! subroutine
     procedure, public :: print    ! subroutine

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
    real(dp), optional :: x
    real(dp), optional :: y
    type(scalar)       :: this
    
    if (present(x))  call this%set_real(x)
    if (present(y))  call this%set_cplx(y)

  end function constructor

  !-------------------------------------------------------------------!
  ! Getter for the real part
  !-------------------------------------------------------------------!
  function get_real(this)

    ! arguments
    class(scalar) :: this
    real(dp)      :: get_real


    ! return value is set
    get_real = this%x

  end function get_real

  !-------------------------------------------------------------------!
  ! getter for the complex part
  !-------------------------------------------------------------------!
  function get_cplx(this)

    ! arguments
    class(scalar) :: this
    real(dp)      :: get_cplx

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

  !-------------------------------------------------------------------!
  ! Print routine to display the contents of the type
  !-------------------------------------------------------------------!
  subroutine print(this)

    ! argumants
    class(scalar) :: this

    write(*,*) this%x, this%y

  end subroutine print

end module scalar_class
