!=====================================================================!
! Module which handles BODIES.
!---------------------------------------------------------------------!
! Has functions to:
! (a) create_body based on supplied state, mass, J, c etc
! (b) update the state variables alone during time-integration
! (d) compute the rotation and angular rate matrices
! (c) 'toString' like function to print the body props (state+attrs)
!---------------------------------------------------------------------!
! Note: 
!---------------------------------------------------------------------!
! Use the interface names to access the functions instead of internal
! functions whose signature may change as per needs.
!=====================================================================!
module bodies

  use types
  use global_constants
  use global_variables
  use utils

  implicit none

  ! Default make everything private
  private 

  ! Expose only the needed variables and functions
  public create_body, set_state 
  public get_rotation
  public get_angrate, get_angrate_dot
  public print_body

  !*******************************************************************!
  ! A common interface for setting (updating the state vectors) during
  ! time-integration. During time-integration only the state of the 
  ! bodies change and it is unnecessary to create new bodies for the 
  ! new state. Therefore the user can call this interface to update the 
  ! state of the supplied  bodies
  !-------------------------------------------------------------------! 
  interface set_state
     module procedure set_body_state
  end interface set_state

  !*******************************************************************!
  ! A common interface for different ways of getting rotation matrix
  !-------------------------------------------------------------------!
  ! (a) get_rotation_from_angles -- > input theta(3)
  ! (b) get_rotation_from_cosines -- > input dir cosines and sines
  ! (c) get_rotation_from_state_vector --> not functional currently
  !*******************************************************************!
  interface get_rotation
     module procedure get_rotation_from_angles_array, &
          &get_rotation_from_angles_vec, get_rotation_from_cosines
     ! get_rotation_from_state_vector
  end interface get_rotation

  !*******************************************************************!
  ! A common interface for different ways of getting rotation matrix
  !-------------------------------------------------------------------!
  ! (a) get_angrate_from_angles -- > input theta(3)
  ! (b) get_angrate_from_cosines -- > input dir cosines and sines
  ! (c) get_angrate_from_state_vector --> not functional currently
  !*******************************************************************!
  interface get_angrate
     module procedure get_angrate_from_angles, get_angrate_from_cosines
     ! get_angrate_from_state_vector
  end interface get_angrate

contains

  !*******************************************************************!
  ! create a body using the supplied state and other parameters
  !*******************************************************************!
  function create_body(mass, re, q, qdot) result(alpha)

    ! inputs
    real(dp), intent(in)        :: mass 
    type(vector), intent(in)    :: re 
    real(dp), intent(in)        :: q(NDOF_PBODY)
    real(dp), intent(in)        :: qdot(NDOF_PBODY)

    ! output
    type(body)                  :: alpha

    ! local variables
    real(dp), dimension(NUM_SPAT_DIM):: theta

    !  sets the state of the body and its time derivatives
    call set_state(q, qdot, alpha)

    ! make a local copy
    theta = q(4:6)

    ! rotation and rate matrices
    alpha%C_mat     = get_rotation(theta)
    alpha%S         = get_angrate(theta)
    alpha%S_dot     = get_angrate_dot(theta, qdot(4:6))

    ! mass of the body
    alpha%mass  = mass

    ! mass moment of intertia (second moment)
    alpha%J = -mass*skew(re)*skew(re)

    ! first moment of mass
    alpha%c  = mass*re

    ! reaction force 
    alpha%fr = zeroV !mass*GRAV !? check

    ! reaction torque
    alpha%gr = zeroV !alpha%J*alpha%omega_dot !? check
    
    ! translational + rotational KE = 0.5mV^2 + 0.5*wJw
    ! (*) between two vectors does a dot product (see utils.f90)
    ! Assuming the body axis to be at the centre of mass of the body
    ! Inertial or body frame?
    alpha%KE = 0.5_dp*(mass * alpha%v *  alpha%v &
         & + alpha%omega*alpha%J* alpha%omega)

    ! potential energy due to position
    ! include strain energy later
    ! Inertial or body frame?
    alpha%PE = mass*GRAV*alpha%r
    
    call print_body(alpha)
    
  end function create_body

  !*******************************************************************!
  ! Takes the body (alpha) and updates/sets the state variables
  !*******************************************************************!
  subroutine set_body_state(q, qdot, alpha)

    real(dp), intent(in)          :: q(NDOF_PBODY)
    real(dp), intent(in)          :: qdot(NDOF_PBODY)
    type(body),intent(inout)      :: alpha

    ! set the state into the body
    alpha%r         = vector(q(1:3))
    alpha%theta     = vector(q(4:6))
    alpha%v         = vector(q(7:9))
    alpha%omega     = vector(q(10:12))

    ! set the time derivatives of state into the body
    alpha%r_dot     = vector(qdot(1:3))
    alpha%theta_dot = vector(qdot(4:6))
    alpha%v_dot     = vector(qdot(7:9))
    alpha%omega_dot = vector(qdot(10:12))
    
    if (ELASTIC) then
       ! set the elastic state variables in the body
!       alpha%qs        = vector(q(13:15))    ! IS_QS:IE_QE
!       alpha%qs_dot    = vector(qdot(13:15)) ! IS_QS:IE_QE
    end if

  end subroutine set_body_state

  !*******************************************************************!
  ! Returns the rotation matrix based on the euler angles
  ! Compute the 3-2-1 Euler angle rotation

  ! C = C1(theta_1)*C2(theta_2)*C3(theta_3)
  !
  ! Input: theta as an array
  ! Output: CMAT of type MATRIX
  !*******************************************************************!
  function get_rotation_from_angles_array(theta) result(CMAT)

    real(dp), intent(in)       :: theta(NUM_SPAT_DIM)    
    type(matrix)               :: CMAT
    real(dp)                   :: c1, c2, c3, s1, s2, s3

    ! Compute the sin/cos of the Euler angles
    c1 = cos(theta(1))
    s1 = sin(theta(1))

    c2 = cos(theta(2))
    s2 = sin(theta(2))

    c3 = cos(theta(3))
    s3 = sin(theta(3))

    CMAT = get_rotation_from_cosines(c1, c2, c3, s1, s2, s3)

  end function get_rotation_from_angles_array
  
  !*******************************************************************!
  ! Returns the rotation matrix based on the euler angles
  ! Compute the 3-2-1 Euler angle rotation
  !
  ! CMAT = C1(theta_1)*C2(theta_2)*C3(theta_3)
  !
  ! Input: theta of type VECTOR
  ! Output: CMAT of type MATRIX
  !*******************************************************************!
  function get_rotation_from_angles_vec(thetain) result(CMAT)

    type(vector), intent(in)   :: thetain
    real(dp)                   :: theta(NUM_SPAT_DIM)
    type(matrix)               :: CMAT
    real(dp)                   :: c1, c2, c3, s1, s2, s3

    ! covert to array form
    theta = array(thetain)

    ! call the method that takes angles array
    CMAT  =  get_rotation_from_angles_array(theta)

  end function get_rotation_from_angles_vec

  !*******************************************************************!
  ! Returns the rotation matrix (euler angles) based on the sines and 
  ! cosines of the euler angles
  !*******************************************************************!
  function get_rotation_from_cosines(c1,c2,c3,s1,s2,s3) result(CMAT)

    real(dp), intent(in)       :: c1, c2, c3, s1, s2, s3
    type(matrix)               :: CMAT

    CMAT = matrix((/ c2*c3, c2*s3, -s2,&
         & s1*s2*c3 - c1*s3, s1*s2*s3 + c1*c3, s1*c2,&
         & c1*s2*c3 + s1*s3, c1*s2*s3 - s1*c3, c1*c2 /))

  end function get_rotation_from_cosines

  ! ********************************************************
  ! routine that returns the ang rate matrix (euler angles)
  ! ********************************************************
  function get_angrate_from_angles(theta) result(SMAT)

    real(dp) :: theta(NUM_SPAT_DIM)    
    type(matrix) :: SMAT
    real(dp)                   :: c1, c2, c3, s1, s2, s3

    ! Compute the sin/cos of the Euler angles
    c1 = cos(theta(1))
    s1 = sin(theta(1))

    c2 = cos(theta(2))
    s2 = sin(theta(2))

    c3 = cos(theta(3))
    s3 = sin(theta(3))

    SMAT = get_angrate_from_cosines(c1, c2, c3, s1, s2, s3)

  end function get_angrate_from_angles

  ! ********************************************************
  ! routine that returns the rotation matrix (euler angles)
  ! based on the sines and cosines of the euler angles
  ! ********************************************************
  function get_angrate_from_cosines(c1,c2,c3,s1,s2,s3) result(SMAT)

    real(dp), intent(in)       :: c1, c2, c3, s1, s2, s3
    type(matrix)               :: SMAT

    SMAT = matrix((/ 1.0_dp, 0.0_dp, -s2, &
         & 0.0_dp,  c1,  s1*c2, &
         & 0.0_dp,  -s1,  c1*c2 /))

  end function get_angrate_from_cosines

  !-----------------------------------------------------------
  ! we use the 3-2-1 Euler angles, the S matrix does not depend
  ! on c3/s3, however we keep these as inputs in case we ever  
  ! want to change the Euler sequence.
  !-----------------------------------------------------------
  function get_angrate_dot(theta, dtheta) result(SMAT_DOT)

    real(dp)     :: theta(NUM_SPAT_DIM), dtheta(NUM_SPAT_DIM)
    type(matrix) :: SMAT_DOT
    real(dp)     :: c1, c2, c3, s1, s2, s3

    ! Compute the sin/cos of the Euler angles
    c1 = cos(theta(1))
    s1 = sin(theta(1))

    c2 = cos(theta(2))
    s2 = sin(theta(2))

    c3 = cos(theta(3))
    s3 = sin(theta(3))
    
    SMAT_DOT= matrix( (/ 0.0_dp,  0.0_dp,  -c2*dtheta(2), &
         & 0.0_dp, -s1*dtheta(1), c1*c2*dtheta(1)-s1*s2*dtheta(2),&
         & 0.0_dp, -c1*dtheta(1), -s1*c2*dtheta(1)-c1*s2*dtheta(2) /))

  end function get_angrate_dot


  ! ********************************************************
  ! routine that prints the state and properties of the body
  ! ********************************************************
  subroutine print_body(alpha)

    use dispmodule

    type(body):: alpha

    call disp('-----------------------------------------------------')
    call disp('---------------------BODY----------------------------')
    call disp('-----------------------------------------------------')
    call disp('')
    call disp('   r          =   ', get_array(alpha%r), SEP=', ', &
         &ORIENT = 'ROW')
    call disp('   theta      =   ', get_array(alpha%theta), SEP=', ',&
         & ORIENT = 'ROW')
    call disp('   v          =   ', get_array(alpha%r), SEP=', ', &
         &ORIENT = 'ROW')
    call disp('   omega      =   ', get_array(alpha%theta), SEP=', ',&
         & ORIENT = 'ROW')
    call DISP('')
    call disp('   r_dot      =   ', get_array(alpha%r), SEP=', ',&
         & ORIENT = 'ROW')
    call disp('   theta_dot  =   ', get_array(alpha%theta), SEP=', ',&
         & ORIENT = 'ROW')
    call disp('   v_dot      =   ', get_array(alpha%r), SEP=', ',&
         & ORIENT = 'ROW')
    call disp('   omega_dot  =   ', get_array(alpha%theta), SEP=', ',&
         &ORIENT = 'ROW')
    call DISP('')
    call disp('   c          =   ', get_array(alpha%c), SEP=', ',&
         & ORIENT = 'ROW')
    call DISP('')
    call disp('   J          =   ', get_matrix(alpha%J))
    call DISP('')
    call disp('   C          =   ', get_matrix(alpha%C_mat))
    call DISP('')
    call disp('   S          =   ', get_matrix(alpha%S))
    call DISP('')
    call disp('   S_dot      =   ', get_matrix(alpha%S_dot))
    call DISP('')
    call disp('   fr         =   ', get_array(alpha%fr), SEP=', ', &
         &ORIENT = 'ROW')

    call disp('   gr         =   ', get_array(alpha%gr), SEP=', ', &
         &ORIENT = 'ROW')
    call DISP('-------------- ENERGY BALANCE ------------------------')
    call disp('   PE         =   ', alpha%PE)
    call disp('   KE         =   ', alpha%KE)
    call disp('   TE         =   ', alpha%KE + alpha%PE)
    call disp('-----------------------------------------------------')
    call disp('')

  end subroutine print_body

end module bodies

!!$ ! ********************************************************
!!$  ! routine that returns the rotation matrix (euler angles)
!!$  ! based on the angles

  ! compiler threw ambiguity 
  ! error so I have commented out this impl for now and will explore
  ! later
!!$  ! ********************************************************
!!$  function get_rotation_from_state_vector(qr) result(CMAT)
!!$
!!$    real(dp), intent(in)       :: qr(NUM_STATES_PER_BODY)    
!!$    type(matrix), intent(out)  :: CMAT
!!$
!!$    CMAT = get_rotation_from_angles(qr(4:6)) ! angles are stored in indices 4:6
!!$
!!$  end function get_rotation_from_state_vector
!!$
!!$ ! ******************************************************************!
!!$  ! Returns the time derivative of rotation matrix (euler angles)
!!$  ! ******************************************************************!
!!$  function get_rotation_dot(theta)
!!$
!!$    real(dp) :: theta(NUM_SPAT_DIM)
!!$    type(matrix) :: get_rotation_dot
!!$
!!$    stop "dummy impl"
!!$
!!$  end function get_rotation_dot
