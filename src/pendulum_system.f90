!=====================================================================!
! Module that sets up pendulum system
!=====================================================================!

module pendulum_class

  use system_class
  use body_class, only  : body
  use joint_class, only : joint

  ! module dependencies
  implicit none

  private
  public :: pendulum
  
  type, extends(system) :: pendulum

     private

     class(body) , allocatable :: bodies(:)
     class(joint), allocatable :: joints(:)

   contains

     procedure :: add_body => add_body
     procedure :: add_joint => add_joint

  end type pendulum
  
  ! interface for the constructor
  interface pendulum
     procedure constructor
  end interface pendulum

contains
  
  !*******************************************************************!
  ! Constructor for pendulum system
  !*******************************************************************!
  function constructor(nbody, njoint) result(this)

     integer        :: nbody
     integer        :: njoint
     type(pendulum) :: this
     
     call this % set_nbody(nbody)
     call this % set_njoint(njoint)

  end function constructor
  
  !*******************************************************************!
  ! Routine to add a body to the system
  !*******************************************************************!
  subroutine add_body(this, bnum, bdy)

    ! arguments
    class(pendulum) :: this
    integer         :: bnum
    class(body)     :: bdy

    print *, 'Added body to the pendulum!'

    allocate(this%bodies(bnum), source = bdy)
    
  end subroutine add_body

  
  !*******************************************************************!
  ! Routine to add a joint to the system
  !*******************************************************************!
  subroutine add_joint(this, jnum, jnt)

    ! arguments
    class(pendulum) :: this
    integer         :: jnum
    class(joint)    :: jnt

    print *, 'Added joint to the pendulum!'

    allocate(this%joints(jnum), source = jnt)
    
  end subroutine add_joint

end module pendulum_class
