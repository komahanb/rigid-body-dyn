!=====================================================================!
! Abstract class that defines joints
!---------------------------------------------------------------------!
! Has functions to:
! (a) getters and setters to joint properties
! (b) abstract constructor to create joints that one should implement 
!     for each joint type downstream extending this class
!---------------------------------------------------------------------!
! Note: 
!---------------------------------------------------------------------!
! Should be able to accomodate both rigid and elatic bodies as we use
! the abstract class for bodies here instead of their respective impls
!=====================================================================!
module joint_class

  use body_class, only: body

  implicit none

  ! Default make everything private
  private 

  ! Expose only the needed variables and functions
  public joint

  ! we create a new datatype called joint
  type, abstract :: joint
     
     private
     
     ! joint number
     integer :: joint_num

     ! spherical, revolute, prismatic, planar
     character(len=10) :: joint_type      
     
     ! the two interacting bodies
     class(body) :: first_body, second_body
     
   contains
     
     procedure :: get_joint_num, set_joint_num
     procedure :: get_joint_type, set_joint_type
     procedure :: get_first_body, set_fist_body
     procedure :: get_second_body, set_second_body

     procedure(create_joint_interface), deferred :: create_joint

  end type joint

  ! inteface for adding joints
  abstract interface

     subroutine create_joint_interface(this, jnum, jtype, &
          & first_body, second_body)

       use body_class, only: body
       import joint

       class(joint) :: this
       integer      :: jnum
       character(len=10) :: jtype
       class(body)  :: first_body, second_body

     end subroutine create_joint_interface

  end interface

contains  

  !*******************************************************************!
  ! Getter for the current joint number
  !*******************************************************************!

  function get_joint_num(this)

    class(joint) :: this
    integer       :: get_joint_num

    get_joint_num = this % joint_num

  end function get_joint_num

  !*******************************************************************!
  ! Setter for the curernt joint number
  !*******************************************************************!

  subroutine set_joint_num(this, jnum)

    class(joint) :: this
    integer       :: jnum

    this % joint_num =  jnum

  end subroutine set_joint_num

  !*******************************************************************!
  ! Getter for the current joint type
  !*******************************************************************!

  function get_joint_type(this)

    class(joint)      :: this
    character(len=10) :: get_joint_type

    get_joint_type = this % joint_type

  end function get_joint_type

  !*******************************************************************!
  ! Setter for the curernt joint type
  !*******************************************************************!

  subroutine set_joint_type(this, jtype)

    class(joint)      :: this
    character(len=10) :: jtype

    this % joint_type =  jtype

  end subroutine set_joint_type
  
  !*******************************************************************!
  ! Getter for the first_body
  !*******************************************************************!
  
  function get_first_body(this)

    class(joint) :: this
    class(body)  :: get_first_body

    get_first_body = this % first_body

  end function get_first_body

  !*******************************************************************!
  ! Setter for the curernt joint type
  !*******************************************************************!

  subroutine set_first_body(this, first_body)

    class(joint):: this
    class(body) :: first_body

    this % first_body =  first_body

  end subroutine set_first_body


  !*******************************************************************!
  ! Getter for the second_body
  !*******************************************************************!
  
  function get_second_body(this)

    class(joint) :: this
    class(body)  :: get_second_body

    get_second_body = this % second_body

  end function get_second_body

  !*******************************************************************!
  ! Setter for the curernt joint type
  !*******************************************************************!

  subroutine set_second_body(this, second_body)

    class(joint):: this
    class(body) :: second_body

    this % second_body =  second_body

  end subroutine set_second_body

end module joint_class

