!=====================================================================================!
! Implementation of the joint_abstract class for speherical joints
! Spherical joint suppots arbitrary rotation about a fixed point
! The effect of elastic deformations are ignored for simplicity
!=====================================================================================!

module spherical_joint_class

  use joint_class
  use body_class
  use types, only: vector

  implicit none

  ! Default make everything private
  private 

  ! Expose only the needed variables and functions
  public spherical_joint

  ! new type for spherical joint
  type, extends(joint) :: spherical_joint

     contains

!     procedure :: joint_residual => joint_residual !spherical_joint_residual

  end type spherical_joint
  
  ! interface for the constructor
  interface spherical_joint
     procedure constructor
  end interface spherical_joint

contains
    
  !*******************************************************************!
  ! Constructor for spherical joint
  !*******************************************************************!
  
  function constructor(jnum, jtype, first_body, second_body) &
       & result(this)

    integer               :: jnum
    character(len=5)      :: jtype
    class(body)           :: first_body, second_body
    type(spherical_joint) :: this

    call this % set_joint_num(jnum)
    call this % set_joint_type(jtype)
    call this % set_first_body(first_body)
    call this % set_second_body(second_body)

  end function constructor

!!$  !*******************************************************************!
!!$  ! Function that implements the joint equations
!!$  !*******************************************************************!
!!$  subroutine spherical_joint_residual(this, residual)
!!$    
!!$    use types, only: vector
!!$    
!!$    type(spherical_joint) :: this
!!$    type(vector) :: residual(2) ! 6 equations (2 in vector form)
!!$    
!!$!    residual(1) = 
!!$
!!$  end subroutine spherical_joint_residual
!!$

end module spherical_joint_class
