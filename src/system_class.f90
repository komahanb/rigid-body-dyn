! Module that defines the system 
module system_class

  ! module dependencies
  use body_class, only  : body
  use joint_class, only : joint

  implicit none

  private
  public :: system

  type system
     
     private

     integer   :: nbody
     integer   :: njoint

   contains

     procedure :: get_nbody
     procedure :: set_nbody

     procedure :: get_njoint
     procedure :: set_njoint

     procedure :: add_body
     procedure :: add_joint

  end type system

!!$  
!!$  ! inteface for creating bodies
!!$  interface add_body
!!$     module procedure add_body_to_system
!!$  end interface add_body
!!$
!!$  interface add_joint
!!$     module procedure add_joint_to_system
!!$  end interface add_joint

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

  !*******************************************************************!
  ! Add body into the system
  !*******************************************************************!

  subroutine add_body(this, bnum, bdy)

    class(system) :: this
    integer       :: bnum
    class(body)   :: bdy

  end subroutine add_body
  
  !*******************************************************************!
  ! Add joint into the system
  !*******************************************************************!
  
  subroutine add_joint(this, jnum, jnt)

    class(system) :: this
    integer       :: jnum
    class(joint)  :: jnt

  end subroutine add_joint

end module system_class
