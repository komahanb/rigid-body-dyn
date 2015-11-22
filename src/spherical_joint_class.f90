!=====================================================================================!
! Implementation of the joint_abstract class for speherical joints
! Spherical joint suppots arbitrary rotation about a fixed point
! The effect of elastic deformations are ignored for simplicity
!=====================================================================================!

module spherical_joint_class

  use global_constants, only : dp, NDOF_PJOINT, NUM_JOINT_EQN
  use joint_class, only : joint
  use body_class, only : body
  use utils, only: operator(*), operator(+), operator(-),&
       & matrix, array, skew
  use types, only: vector

  implicit none

  ! Default make everything private
  private 

  ! Expose only the needed variables and functions
  public spherical_joint
  
  ! new type for spherical joint
  type, extends(joint) :: spherical_joint

   contains

     procedure :: get_residual => get_joint_residual

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

  !*******************************************************************!
  ! Function that implements the joint equations
  !*******************************************************************!

  function get_joint_residual(this) result (residual)
    
    class(spherical_joint) :: this
    
    type(vector) :: res_vec(NUM_JOINT_EQN) ! 6 equations (2 in vector form)
    real(dp)     :: residual(NDOF_PJOINT) ! scalar form
    
    !res_vec(1) = 
    !res_vec(2) = 
    
    !-----------------------------------------------------------------!
    ! convert vector to array form
    !-----------------------------------------------------------------!
    
    residual = array( res_vec )
    
  end function get_joint_residual


end module spherical_joint_class
