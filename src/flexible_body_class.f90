module flexible_body_class

  ! module references
  use global_constants
  use global_variables
  use types
  use utils
  use body_class
  use rigid_body_class
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
  interface flexible_body

     procedure constructor ! add constructor

  end interface flexible_body

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
    real(dp), optional, intent(in) :: mass 
    real(dp), optional, intent(in) :: c(NUM_SPAT_DIM)
    real(dp), optional, intent(in) :: J(NUM_SPAT_DIM, NUM_SPAT_DIM)
    real(dp), optional, intent(in) :: fr(NUM_SPAT_DIM)
    real(dp), optional, intent(in) :: gr(NUM_SPAT_DIM)

    real(dp), optional, intent(in) :: q(NDOF_PBODY)
    real(dp), optional, intent(in) :: qdot(NDOF_PBODY)

    real(dp), optional, intent(in) :: qs(NUM_SPAT_DIM)
    real(dp), optional, intent(in) :: qs_dot(NUM_SPAT_DIM)
    real(dp), optional, intent(in) :: qs_double_dot(NUM_SPAT_DIM)

    real(dp), optional, intent(in) :: f(NUM_SPAT_DIM)

    real(dp), optional, intent(in) :: p(NUM_SPAT_DIM, NUM_SPAT_DIM)
    real(dp), optional, intent(in) :: h(NUM_SPAT_DIM, NUM_SPAT_DIM)
    real(dp), optional, intent(in) :: K(NUM_SPAT_DIM, NUM_SPAT_DIM)
    real(dp), optional, intent(in) :: M(NUM_SPAT_DIM, NUM_SPAT_DIM)

    ! input/output
    type(flexible_body):: this
    !-----------------------------------------------------------------!
    ! Inertial properties of the body
    !-----------------------------------------------------------------!

    ! mass of the body
    if (present(mass)) this%mass = mass

    ! moment of inertia in body-fixed frame
    if (present(J)) this%J = matrix(J)       !-mass*skew(re)*skew(re)

    !first moment of inertia in body-fixed frame: mass*(cg location)
    if (present(c)) this%c = vector(c)       ! mass*re

    !-----------------------------------------------------------------!
    ! set the state into the body
    !-----------------------------------------------------------------!
    if (present(q)) then
       this%r         = vector(q(1:3))
       this%theta     = vector(q(4:6))
       this%v         = vector(q(7:9))
       this%omega     = vector(q(10:12))
    end if

    !-----------------------------------------------------------------!
    ! set the time derivatives of state into the body
    !-----------------------------------------------------------------!

    if (present(qdot)) then
       this%r_dot     = vector(qdot(1:3))
       this%theta_dot = vector(qdot(4:6))
       this%v_dot     = vector(qdot(7:9))
       this%omega_dot = vector(qdot(10:12))
    end if

    !-----------------------------------------------------------------!
    ! update the rotation and angular rate matrices
    !-----------------------------------------------------------------!

    this%C_mat     = get_rotation(this%theta)
    this%S         = get_angrate(this%theta)
    this%S_dot     = get_angrate_dot(this%theta, this%theta_dot)

    !-----------------------------------------------------------------!
    ! update the new direction of the gravity vector in body frame
    !-----------------------------------------------------------------!

    this%g = this%C_mat*GRAV

    !-----------------------------------------------------------------!
    ! Mechanical Energy 
    !-----------------------------------------------------------------!

    ! update the kinetic energy
    this%KE = 0.5_dp*(this%mass * this%v *  this%v &
         & + this%omega*this%J* this%omega) ! + coupling term

    ! update potential energy
    this%PE = this%mass*this%g*this%r

    !-----------------------------------------------------------------!
    ! Joint reactions
    !-----------------------------------------------------------------!

    ! reaction force 
    if (present(fr)) this%fr    = vector(fr)   

    ! reaction torque
    if (present(gr)) this%gr    = vector(gr)

    if (this%mass .eq. ZERO) stop "ERROR: Body with zero mass!"


    if (present(qs))  this%qs = vector(qs)
    if (present(qs_dot))  this%qs_dot = vector(qs_dot)
    if (present(qs_double_dot))  this%qs_double_dot = &
         &vector(qs_double_dot)

    if (present(f))  this%f = vector(f)

    if (present(p))  this%p = matrix(p)
    if (present(h))  this%h = matrix(h)
    if (present(K))  this%K = matrix(K)
    if (present(M))  this%M = matrix(M)

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
