!=====================================================================!
! Abstract class that defines joints
!---------------------------------------------------------------------!
! Has functions to:
!
! (a) getters and setters to joint properties
! (b) abstract constructor to create joints that one should implement 
!     for each joint type downstream extending this class
!---------------------------------------------------------------------!
! Note: 
!---------------------------------------------------------------------!
!
! Should be able to accomodate both rigid and elatic bodies as we use
! the abstract class for bodies here instead of their respective impls
!
! We formulate six equations for each type of joint irrespective of
! whether or not some are trivial. This is because:
!
! (1) the structure of the equations is not affected by the joints
! (2) helps exploit sparsity
!
!=====================================================================!
module joint_class

  use body_class, only: body
  use types, only: vector

  implicit none

  ! Default make everything private
  private 

  ! Expose only the needed variables and functions
  public joint

  ! we create a new datatype called joint
  type, abstract :: joint

     private

     integer :: joint_num           ! joint number
     character(len=5) :: joint_type ! spherical, revolute, prismatic, planar
     ! could just use ra and c_a but will leave it for now
     class(body), allocatable :: first_body, second_body ! the two interacting bodies

   contains

     procedure :: get_joint_num, set_joint_num
     procedure :: get_joint_type, set_joint_type
     procedure :: get_first_body , set_first_body
     procedure :: get_second_body, set_second_body

     procedure(get_residual_interface), deferred :: get_residual

  end type joint


  ! an abstract-interface for finding residual of the body
  abstract interface

     function get_residual_interface(this) result(res)
       use global_constants, only : dp, NDOF_PJOINT
       import joint
       class(joint) :: this
       real(dp)    :: res(NDOF_PJOINT)
     end function get_residual_interface
     
  end interface
  
contains  

  !*******************************************************************!
  ! Getter for the current joint number
  !*******************************************************************!

  function get_joint_num(this)

    class(joint)  :: this
    integer       :: get_joint_num

    get_joint_num = this % joint_num

  end function get_joint_num

  !*******************************************************************!
  ! Setter for the curernt joint number
  !*******************************************************************!

  subroutine set_joint_num(this, jnum)

    class(joint)  :: this
    integer       :: jnum

    this % joint_num =  jnum

  end subroutine set_joint_num

  !*******************************************************************!
  ! Getter for the current joint type
  !*******************************************************************!

  function get_joint_type(this)

    class(joint)     :: this
    character(len=5) :: get_joint_type

    get_joint_type = this % joint_type

  end function get_joint_type

  !*******************************************************************!
  ! Setter for the curernt joint type
  !*******************************************************************!

  subroutine set_joint_type(this, jtype)

    class(joint)     :: this
    character(len=5) :: jtype

    this % joint_type =  jtype

  end subroutine set_joint_type

  !*******************************************************************!
  ! Getter for the first_body
  !*******************************************************************!

  function get_first_body(this)

    class(joint) :: this
    class(body), allocatable  :: get_first_body

    allocate(get_first_body, source = this%first_body)

  end function get_first_body

  !*******************************************************************!
  ! Setter for the curernt joint type
  !*******************************************************************!

  subroutine set_first_body(this, first_body)

    class(joint):: this
    class(body) :: first_body

    allocate(this%first_body, source = first_body)

  end subroutine set_first_body


  !*******************************************************************!
  ! Getter for the second_body
  !*******************************************************************!

  function get_second_body(this)

    class(joint) :: this
    class(body), allocatable :: get_second_body

    allocate(get_second_body, source = this%second_body)

  end function get_second_body

  !*******************************************************************!
  ! Setter for the curernt joint type
  !*******************************************************************!

  subroutine set_second_body(this, second_body)

    class(joint):: this
    class(body) :: second_body

    allocate(this % second_body, source=second_body)

  end subroutine set_second_body

end module joint_class
