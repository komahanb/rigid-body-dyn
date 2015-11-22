!=====================================================================!
! Implementation of the joint_abstract class for speherical joints
! Spherical joint suppots arbitrary rotation about a fixed point
! The effect of elastic deformations are ignored for simplicity
!=====================================================================!

module spherical_joint_class

  use global_constants, only : dp, NDOF_PJOINT, NUM_JOINT_EQN
  use global_variables
  use joint_class, only : joint
  use body_class, only : body
  use rigid_body_class
  use utils, only: operator(*), operator(+), operator(-),&
       & matrix, array, skew, trans
  use types, only: vector, matrix

  implicit none

  ! Default make everything private
  private 

  ! Expose only the needed variables and functions
  public spherical_joint
  
  ! new type for spherical joint
  type, extends(joint) :: spherical_joint

   contains

     ! implementations for deferred functions and routines from abstract class

     procedure :: get_residual => get_joint_residual
     procedure :: get_jacobian => get_joint_jacobian

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
    
    real(dp)     :: residual(NDOF_PJOINT) ! scalar form
    type(vector) :: res_vec(NUM_JOINT_EQN) ! 6 equations (2 in vector form)

    ! local variables
    type(vector) :: ra1, ra, rb, rb1
    type(matrix) :: C_a, C_b

    class(body), allocatable   :: bodyA, bodyB
    
    ! set the local variables

    allocate(bodyA, source = this % get_first_body() )
    allocate(bodyB, source = this % get_second_body() )
    
    ra = bodyA % get_r()
    rb = bodyB % get_r()

    ra1 = bodyA % get_joint_location()
    rb1 = bodyB % get_joint_location()

    C_a = bodyA % get_rotation()
    C_b = bodyB % get_rotation()

    ! Joint equations
    
    res_vec(1)  = ra + trans(C_a) * ra1 - rb - trans(C_b) * rb1 !force
    res_vec(2)  = zeroV ! torque

    !-----------------------------------------------------------------!
    ! convert vector to array form
    !-----------------------------------------------------------------!
    
    residual = array(res_vec)
    
  end function get_joint_residual

  
  !*******************************************************************!
  ! Function that implements the jacobian for the joint equations
  !*******************************************************************!
  
  function get_joint_jacobian(this) result (jacobian)

    class(spherical_joint) :: this
    real(dp)     :: jacobian(NDOF_PJOINT, NDOF_PJOINT)      ! 6 x 6 scalar form

    type(matrix) :: jac_block(NUM_JOINT_EQN, NUM_JOINT_EQN) ! 2 X 2 in matrix-block form

!!$    ! local variables
!!$    type(vector) :: ra1, ra, rb, rb1
!!$    type(matrix) :: C_a, C_b
!!$
!!$    class(body), allocatable   :: bodyA, bodyB

    ! set the local variables

!!$    allocate(bodyA, source = this % get_first_body() )
!!$    allocate(bodyB, source = this % get_second_body() )
!!$    
!!$    ra = bodyA % get_r()
!!$    rb = bodyB % get_r()
!!$
!!$    ra1 = bodyA % get_joint_location()
!!$    rb1 = bodyB % get_joint_location()
!!$
!!$    C_a = bodyA % get_rotation()
!!$    C_b = bodyB % get_rotation()

    ! Joint jacobian

!!$    jac_block(1,1) = 
!!$    jac_block(2,1) = 
!!$    jac_block(1,2) = 
!!$    jac_block(2,2) = 

    !-----------------------------------------------------------------!
    ! convert block to primitive matrix form
    !-----------------------------------------------------------------!

    jacobian = matrix(jac_block)

  end function get_joint_jacobian

end module spherical_joint_class
