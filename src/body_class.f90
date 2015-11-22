!=====================================================================!
! An abstract class to handle bodies
!---------------------------------------------------------------------!
! This class should be extended to implement
! (a) rigid body
! (b) flexible body
!---------------------------------------------------------------------!
! This abstract class has:
!
! Getters and setters that are relevant to the variables in the type
! Method to initialize the body
!=====================================================================!

module body_class
  
  ! module dependencies
  use global_constants, only : dp

  implicit none

  ! Default make everything private
  private 

  ! Expose only the needed variables and functions
  public :: body

  ! Abtract that that subclass must extend
  type, abstract :: body

     integer          :: body_num     ! body number
     character(len=5) :: body_type    ! rod, bar, sphere, plate

   contains 

     procedure :: get_body_type, set_body_type     
     procedure :: get_body_num, set_body_num

     procedure(get_residual_interface), deferred :: get_residual
     procedure(print_body_interface), deferred :: print

  end type body
  
  ! an abstract-interface for finding residual of the body
  abstract interface

     function get_residual_interface(this) result(res)
       use global_constants, only : dp, NDOF_PBODY
       import body
       class(body) :: this
       real(dp)    :: res(NDOF_PBODY)
     end function get_residual_interface
     
  end interface

  ! an abstract-interface for print implementations
  abstract interface

     subroutine print_body_interface(this)
       import body
       class(body) :: this
     end subroutine print_body_interface

  end interface



contains

  !*******************************************************************!
  ! Getter for the current body number
  !*******************************************************************!

  function get_body_num(this)

    class(body) :: this
    integer     :: get_body_num

    get_body_num = this % body_num

  end function get_body_num

  !*******************************************************************!
  ! Setter for the curernt body number
  !*******************************************************************!

  subroutine set_body_num(this, bnum)

    class(body) :: this
    integer     :: bnum

    this % body_num =  bnum

  end subroutine set_body_num

  !*******************************************************************!
  ! Getter for the current body type
  !*******************************************************************!

  function get_body_type(this)

    class(body)      :: this
    character(len=5) :: get_body_type

    get_body_type = this % body_type

  end function get_body_type

  !*******************************************************************!
  ! Setter for the curernt body type
  !*******************************************************************!

  subroutine set_body_type(this, btype)

    class(body)      :: this
    character(len=5) :: btype

    this % body_type =  btype

  end subroutine set_body_type

end module body_class
