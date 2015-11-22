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

     procedure(get_r_interface), deferred        :: get_r
     procedure(get_rotation_interface), deferred :: get_rotation
     procedure(get_joint_location_interface), deferred :: &
          & get_joint_location

     procedure(get_residual_interface), deferred :: get_residual
     procedure(get_jacobian_interface), deferred :: get_jacobian

     procedure(print_body_interface), deferred :: print

  end type body
  
  abstract interface

     ! an abstract-interface for finding residual of the body
     function get_residual_interface(this) result(res)
       use global_constants, only : dp, NDOF_PBODY
       import body
       class(body) :: this
       real(dp)    :: res(NDOF_PBODY)
     end function get_residual_interface
     
     ! an abstract-interface for finding jacobian of the joint
     
     function get_jacobian_interface(this) result(jacobian)
       use global_constants, only : dp, NDOF_PBODY
       import body
       class(body)  :: this
       real(dp)     :: jacobian(NDOF_PBODY, NDOF_PBODY) ! maybe sparse format
     end function get_jacobian_interface

     ! an abstract-interface for print implementations

     subroutine print_body_interface(this)
       import body
       class(body) :: this
     end subroutine print_body_interface

     ! an abstract-interface for implementing rotation matrix
     function get_rotation_interface(this)
       use types, only : matrix
       import body
       class(body) :: this
       type(matrix) :: get_rotation_interface
     end function get_rotation_interface

     ! an abstract-interface for returning the position of body axis
     function get_r_interface(this)
       use types, only : vector
       import body
       class(body)  :: this
       type(vector) :: get_r_interface
     end function get_r_interface

     !inteface for getting the joint location on the body
     function get_joint_location_interface(this)
       use types, only : vector
       import body
       class(body)  :: this
       type(vector) :: get_joint_location_interface
     end function get_joint_location_interface

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
