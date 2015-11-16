! Module that defines the system 
module system_class

  implicit none

  private
  public :: system

  type, abstract :: system

     private

     integer   :: num_body

   contains

     procedure :: get_num_body
     procedure :: set_num_body

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

  function get_num_body(this)

    class(system) :: this
    integer       :: get_num_body

    get_num_body = this % num_body

  end function get_num_body

  !*******************************************************************!
  ! Setter for number of bodies
  !*******************************************************************!

  subroutine set_num_body(this, bnum)

    class(system) :: this
    integer       :: bnum

    this % num_body =  bnum

  end subroutine set_num_body


end module system_class
