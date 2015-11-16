! Module that defines the system 
module system_class

  implicit none

  private
  public :: system

  type, abstract :: system

     private

     integer   :: nbody

   contains

     procedure :: get_nbody
     procedure :: set_nbody

     ! deferred procedure
     procedure(add_body_interface), deferred :: add_body

  end type system

  ! inteface for creating bodies
  abstract interface

     subroutine add_body_interface(this, bnum, bdy)

       use body_class, only : body

       import system

       class(system):: this
       integer      :: bnum
       class(body)  :: bdy

     end subroutine add_body_interface

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


end module system_class
