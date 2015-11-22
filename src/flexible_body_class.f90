!=====================================================================!
! Flexible body implementation of the abstract body class extending 
! rigid body class
!---------------------------------------------------------------------!
! Has functions to:
! (a) create flexible body based on the supplied parameters
! (b) update the flexible body state variables during time-marching
! (d) compute the rotation and angular rate matrices of the body
! (c) 'toString' like function to print the body props (state+attrs)
!=====================================================================!

module flexible_body_class

  ! module references
  use rigid_body_class
  use global_constants
  use global_variables
  use types
  use utils
  use body_class
  use dispmodule, only : disp


  ! module options
  implicit none
  private
  public :: flexible_body, print_flexible_body

  ! module definitions
  type, extends(rigid_body) :: flexible_body

     !----------------------------------------------------------------!
     ! elatic state variables and time derivatives
     !----------------------------------------------------------------!
  type(vector) :: qs
  type(vector) :: qs_dot
  type(vector) :: qs_double_dot

  type(matrix) :: p              ! 
  type(matrix) :: h              ! 

  type(matrix) :: K              ! stiffness matrix
  type(matrix) :: M              ! mass matrix

  type(vector) :: f              ! elastic force

end type flexible_body

! interfaces
!  interface flexible_body

!     procedure constructor ! add constructor

!  end interface flexible_body

contains

  !*******************************************************************!
  ! create a body using the supplied parameters
  !-------------------------------------------------------------------!
  ! mass   : mass of the body
  ! c      : first moment of mass 
  ! J      : second moment of mass
  ! fr     : reaction force
  ! gr     : reaction torque
  ! q, qdot: state vector and time derivatives
  ! qddot  : second time derivative of the state (used in elastic only)
  !*******************************************************************!
function constructor(mass, c, J, fr, gr, q, qdot, qs, &
    & qs_dot, qs_double_dot, f, p, h, K, M) result(this)

    ! inputs
    real(dp), intent(in) :: mass 
    real(dp), intent(in) :: c(NUM_SPAT_DIM)
    real(dp), intent(in) :: J(NUM_SPAT_DIM, NUM_SPAT_DIM)
    real(dp), intent(in) :: fr(NUM_SPAT_DIM)
    real(dp), intent(in) :: gr(NUM_SPAT_DIM)

    real(dp), intent(in) :: q(NDOF_PBODY)
    real(dp), intent(in) :: qdot(NDOF_PBODY)

    real(dp), intent(in) :: qs(NUM_SPAT_DIM)
    real(dp), intent(in) :: qs_dot(NUM_SPAT_DIM)
    real(dp), intent(in) :: qs_double_dot(NUM_SPAT_DIM)

    real(dp), intent(in) :: f(NUM_SPAT_DIM)

    real(dp), intent(in) :: p(NUM_SPAT_DIM, NUM_SPAT_DIM)
    real(dp), intent(in) :: h(NUM_SPAT_DIM, NUM_SPAT_DIM)
    real(dp), intent(in) :: K(NUM_SPAT_DIM, NUM_SPAT_DIM)
    real(dp), intent(in) :: M(NUM_SPAT_DIM, NUM_SPAT_DIM)

    ! input/output
    type(flexible_body)  :: this

!!$    !this = rigid_body(mass, c, J, fr, gr, q, qdot)
!!$    
!!$    !-----------------------------------------------------------------!
!!$    ! Inertial properties of the body
!!$    !-----------------------------------------------------------------!
!!$
!!$    ! mass of the body
!!$    this%mass = mass
!!$
!!$    ! moment of inertia in body-fixed frame
!!$    this%J = matrix(J)       !-mass*skew(re)*skew(re)
!!$
!!$    !first moment of inertia in body-fixed frame: mass*(cg location)
!!$    this%c = vector(c)       ! mass*re
!!$
!!$    !-----------------------------------------------------------------!
!!$    ! set the state into the body
!!$    !-----------------------------------------------------------------!
!!$    
!!$       this%r         = vector(q(1:3))
!!$       this%theta     = vector(q(4:6))
!!$       this%v         = vector(q(7:9))
!!$       this%omega     = vector(q(10:12))
!!$    
!!$
!!$    !-----------------------------------------------------------------!
!!$    ! set the time derivatives of state into the body
!!$    !-----------------------------------------------------------------!
!!$
!!$    
!!$       this%r_dot     = vector(qdot(1:3))
!!$       this%theta_dot = vector(qdot(4:6))
!!$       this%v_dot     = vector(qdot(7:9))
!!$       this%omega_dot = vector(qdot(10:12))
!!$    
!!$
!!$    !-----------------------------------------------------------------!
!!$    ! update the rotation and angular rate matrices
!!$    !-----------------------------------------------------------------!
!!$
!!$    this%C_mat     = get_rotation(this%theta)
!!$    this%S         = get_angrate(this%theta)
!!$    this%S_dot     = get_angrate_dot(this%theta, this%theta_dot)
!!$
!!$    !-----------------------------------------------------------------!
!!$    ! update the new direction of the gravity vector in body frame
!!$    !-----------------------------------------------------------------!
!!$
!!$    this%g = this%C_mat*GRAV
!!$
!!$    !-----------------------------------------------------------------!
!!$    ! Mechanical Energy 
!!$    !-----------------------------------------------------------------!
!!$
!!$    ! update the kinetic energy
!!$    this%KE = 0.5_dp*(this%mass * this%v *  this%v &
!!$         & + this%omega*this%J* this%omega) ! + coupling term
!!$
!!$    ! update potential energy
!!$    this%PE = this%mass*this%g*this%r
!!$
!!$    !-----------------------------------------------------------------!
!!$    ! Joint reactions
!!$    !-----------------------------------------------------------------!
!!$
!!$    ! reaction force 
!!$    this%fr    = vector(fr)   
!!$
!!$    ! reaction torque
!!$    this%gr    = vector(gr)
!!$
!!$    if (this%mass .eq. ZERO) stop "ERROR: Body with zero mass!"
!!$
!!$
!!$    this%qs = vector(qs)
!!$    this%qs_dot = vector(qs_dot)
!!$    this%qs_double_dot = vector(qs_double_dot)
!!$
!!$    this%f = vector(f)
!!$
!!$    this%p = matrix(p)
!!$    this%h = matrix(h)
!!$    this%K = matrix(K)
!!$    this%M = matrix(M)


 stop "dummy impl"

end function constructor

!*******************************************************************!
! routine that prints the state and properties of the flexible body
!******************************************************************!
subroutine print_flexible_body(this)

 use dispmodule

 class(flexible_body):: this

 print*,"Implement please!!"

end subroutine print_flexible_body

end module flexible_body_class
