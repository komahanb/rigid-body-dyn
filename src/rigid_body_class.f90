module rigid_body_class

  ! module references
  use body_class
  use global_constants
  use global_variables
  use types
  use utils
  use dispmodule, only : disp

  ! module options
  implicit none
  private
  public :: rigid_body, print_rigid_body

  ! module definitions
  type, extends(body) :: rigid_body

  end type rigid_body

!*******************************************************************!
! A common interface for creating bodies and updating existing ones
!-------------------------------------------------------------------!
! The inputs are described in the methods under the interface. One 
! can use the same interface for updating existing bodies.
! Shortly: use this to create and update bodies
!*******************************************************************!

interface rigid_body

   procedure constructor

end interface rigid_body

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
  function constructor(mass, c, J, fr, gr, q, qdot) result(this)

    ! inputs
    real(dp), optional, intent(in) :: mass 
    real(dp), optional, intent(in) :: c(NUM_SPAT_DIM)
    real(dp), optional, intent(in) :: J(NUM_SPAT_DIM, NUM_SPAT_DIM)
    real(dp), optional, intent(in) :: fr(NUM_SPAT_DIM)
    real(dp), optional, intent(in) :: gr(NUM_SPAT_DIM)

    real(dp), optional, intent(in) :: q(NDOF_PBODY)
    real(dp), optional, intent(in) :: qdot(NDOF_PBODY)

    ! input/output
    type(rigid_body) :: this

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

    !call print_body(this)

  end function constructor
  
  !*******************************************************************!
  ! routine that prints the state and properties of the body
  !******************************************************************!
  subroutine print_rigid_body(this)

    class(rigid_body):: this

    call disp('======================================================')
    call disp('---------------------BODY----------------------------')
    call disp('======================================================')
    call disp('')
    call disp('   > r          =   ', array(this%r), SEP=', ', &
         &ORIENT = 'ROW')
    call disp('   > theta      =   ', array(this%theta), SEP=', ',&
         & ORIENT = 'ROW')
    call disp('   > v          =   ', array(this%v), SEP=', ', &
         &ORIENT = 'ROW')
    call disp('   > omega      =   ', array(this%omega), SEP=', ',&
         & ORIENT = 'ROW')
    call DISP('')
    call disp('   > r_dot      =   ', array(this%r_dot), SEP=', ',&
         & ORIENT = 'ROW')
    call disp('   > theta_dot  =   ', array(this%theta_dot), SEP=', ',&
         & ORIENT = 'ROW')
    call disp('   > v_dot      =   ', array(this%v_dot), SEP=', ',&
         & ORIENT = 'ROW')
    call disp('   > omega_dot  =   ', array(this%omega_dot), SEP=', ',&
         &ORIENT = 'ROW')
    call DISP('')
    call disp('   > c          =   ', array(this%c), SEP=', ',&
         & ORIENT = 'ROW')
    call DISP('')
    call disp('   > J          =   ', matrix(this%J))
    call DISP('')
    call disp('   > Rot Mat    =   ', matrix(this%C_mat))
    call DISP('')
    call disp('   > Angrt Mat  =   ', matrix(this%S))
    call DISP('')
    call disp('   > S_dot      =   ', matrix(this%S_dot))
    call DISP('')
    call disp('   > fr         =   ', array(this%fr), SEP=', ', &
         &ORIENT = 'ROW')

    call disp('   > gr         =   ', array(this%gr), SEP=', ', &
         &ORIENT = 'ROW')
    call disp('')
    call disp('   > gravity    =   ', array(this%g), SEP=', ',&
         & ORIENT = 'ROW')
    call disp('')
    call disp('   > Pot Energy =   ', this%PE)
    call disp('   > Kin Energy =   ', this%KE)
    call disp('   > Tot Energy =   ', this%KE + this%PE)
    call disp('======================================================')
  end subroutine print_rigid_body

end module rigid_body_class
