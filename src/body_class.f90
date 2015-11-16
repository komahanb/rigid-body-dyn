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
!
!=====================================================================!

module body_class

  implicit none

  ! Default make everything private
  private 
  
  ! Expose only the needed variables and functions
  public body
  
  ! abtract that that subclass must extend
  type, abstract :: body
     
     private 

     ! body number
     integer :: body_num

     ! rod, bar, sphere, plate
     character(len=10) :: body_type      

   contains

     procedure :: get_body_num, set_body_num
     procedure :: get_body_type, set_body_type

  end type body

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

    class(body)       :: this
    character(len=10) :: get_body_type

    get_body_type = this % body_type

  end function get_body_type

  !*******************************************************************!
  ! Setter for the curernt body type
  !*******************************************************************!

  subroutine set_body_type(this, btype)

    class(body)       :: this
    character(len=10) :: btype

    this % body_type =  btype

  end subroutine set_body_type

end module body_class
