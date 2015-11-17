!=====================================================================!
! An abstract moduel that defines the system 
!=====================================================================!
module system_class

  ! module dependencies
  implicit none

  private
  public :: system

  type, abstract:: system

     private

     integer   :: nbody
     integer   :: njoint

   contains

     procedure :: get_nbody
     procedure :: set_nbody

     procedure :: get_njoint
     procedure :: set_njoint

     procedure(iadd_body), deferred :: add_body
     procedure(iadd_joint), deferred :: add_joint

  end type system
  
  abstract interface

     !****************************************************************!
     ! Add body into the system
     !****************************************************************!

     subroutine iadd_body(this, bnum, bdy)

       use body_class, only  : body
       import system

       class(system) :: this
       integer       :: bnum
       class(body)   :: bdy

     end subroutine iadd_body

  end interface

  abstract interface 

     !****************************************************************!
     ! Add joint into the system
     !****************************************************************!

     subroutine iadd_joint(this, jnum, jnt)

       use joint_class, only : joint
       import system

       class(system) :: this
       integer       :: jnum
       class(joint)  :: jnt

     end subroutine iadd_joint

  end interface

contains

  !*******************************************************************!
  ! Getter for number of bodies
  !*******************************************************************!

  function get_nbody(this)

    class(system) :: this
    integer       :: get_nbody

    get_nbody = this % nbody

  end function get_nbody

  !*******************************************************************!
  ! Setter for number of bodies
  !*******************************************************************!

  subroutine set_nbody(this, bnum)

    class(system) :: this
    integer       :: bnum

    this % nbody =  bnum

  end subroutine set_nbody

   
  !*******************************************************************!
  ! Getter for number of bodies
  !*******************************************************************!
  
  function get_njoint(this)

    class(system) :: this
    integer       :: get_njoint

    get_njoint = this % njoint

  end function get_njoint

  !*******************************************************************!
  ! Setter for number of bodies
  !*******************************************************************!

  subroutine set_njoint(this, jnum)

    class(system) :: this
    integer       :: jnum

    this % njoint =  jnum

  end subroutine set_njoint

end module system_class
