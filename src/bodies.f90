! module that defines the configuration of a single body
module bodies

  use types
  use global_constants
  use global_variables
  use utils

  implicit none

  ! Default make everything private
  private 

  ! Expose only the needed variables and functions
  public set_state, create_body, print_body
  !  public get_rotation, get_rotation_dot
  !  public get_angular_rate, get_angular_rate_dot

  interface set_state
     module procedure set_body_state
  end interface set_state

  !------------------------------------------------------------------
  ! A common interface for different ways of getting rotation matrix
  ! get_rotation_from_angles -- > user supplies theta(3)
  ! get_rotation_from_cosines -- > use supplies direction cosines and sines
  ! get_rotation_from_state_vector --> compiler threw ambiguity error, so I have commented this impl
  !------------------------------------------------------------------
  interface get_rotation
     module procedure get_rotation_from_angles, get_rotation_from_cosines
     ! get_rotation_from_state_vector
  end interface get_rotation

  !------------------------------------------------------------------
  ! A common interface for different ways of getting angular_rate matrix
  ! get_angular_rate_from_angles -- > user supplies theta(3)
  ! get_angular_rate_from_cosines -- > use supplies direction cosines and sines
  ! get_angular_rate_from_state_vector --> compiler threw ambiguity error, so I have commented this impl
  !------------------------------------------------------------------
  interface get_angular_rate
     module procedure get_angular_rate_from_angles, get_angular_rate_from_cosines
     ! get_angular_rate_from_state_vector
  end interface get_angular_rate

contains

  !****************************************************
  ! Takes the body (alpha) and updates/sets the state variables
  !****************************************************
  subroutine set_body_state(q, qdot, alpha)

    real(dp), intent(in)          :: q(NUM_STATES_PER_BODY)
    real(dp), intent(in)          :: qdot(NUM_STATES_PER_BODY)
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

!!$    if (NUM_ELAST_EQN .gt. 0 .and. elastic) then
!!$
!!$       ! set the elastic state variables in the body
!!$       alpha%qs        = vector(qdot(13:15)) ! IS_QS:IE_QE
!!$       alpha%qs_dot    = vector(qdot(13:15)) ! IS_QS:IE_QE
!!$
!!$    end if

  end subroutine set_body_state
!!$
!!$ ! ********************************************************
!!$  ! routine that returns the rotation matrix (euler angles)
!!$  ! based on the angles
!!$  ! ********************************************************
!!$  function get_rotation_from_state_vector(qr) result(CMAT)
!!$
!!$    real(dp), intent(in)       :: qr(NUM_STATES_PER_BODY)    
!!$    type(matrix), intent(out)  :: CMAT
!!$
!!$    CMAT = get_rotation_from_angles(qr(4:6)) ! angles are stored in indices 4:6
!!$
!!$  end function get_rotation_from_state_vector

  ! ********************************************************
  ! routine that returns the rotation matrix (euler angles)
  ! based on the angles
  ! ********************************************************
  function get_rotation_from_angles(theta) result(CMAT)

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

  end function get_rotation_from_angles

  ! ********************************************************
  ! routine that returns the rotation matrix (euler angles)
  ! based on the sines and cosines of the euler angles
  ! ********************************************************
  function get_rotation_from_cosines(c1, c2, c3, s1, s2, s3) result(CMAT)

    real(dp), intent(in)       :: c1, c2, c3, s1, s2, s3
    type(matrix)               :: CMAT

    CMAT = matrix((/ c2*c3, c2*s3, -s2,&
         & s1*s2*c3 - c1*s3, s1*s2*s3 + c1*c3, s1*c2,&
         & c1*s2*c3 + s1*s3, c1*s2*s3 - s1*c3, c1*c2 /))

  end function get_rotation_from_cosines


  ! **************************************************************************
  ! routine that returns the time derivative of rotation matrix (euler angles)
  ! *************************************************************************
  function get_rotation_dot(theta)

    real(dp) :: theta(NUM_SPAT_DIM)
    type(matrix) :: get_rotation_dot

    stop "dummy impl"

  end function get_rotation_dot


  ! ********************************************************
  ! routine that returns the ang rate matrix (euler angles)
  ! ********************************************************
  function get_angular_rate_from_angles(theta) result(SMAT)

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

    SMAT = get_angular_rate_from_cosines(c1, c2, c3, s1, s2, s3)

  end function get_angular_rate_from_angles

  ! ********************************************************
  ! routine that returns the rotation matrix (euler angles)
  ! based on the sines and cosines of the euler angles
  ! ********************************************************
  function get_angular_rate_from_cosines(c1, c2, c3, s1, s2, s3) result(SMAT)

    real(dp), intent(in)       :: c1, c2, c3, s1, s2, s3
    type(matrix)               :: SMAT

    SMAT = matrix((/ 1.0_dp, 0.0_dp, -s2, &
         & 0.0_dp,  c1,  s1*c2, &
         & 0.0_dp,  -s1,  c1*c2 /))

  end function get_angular_rate_from_cosines

  !-----------------------------------------------------------
  ! we use the 3-2-1 Euler angles, the S matrix does not depend
  ! on c3/s3, however we keep these as inputs in case we ever  
  ! want to change the Euler sequence.
  !-----------------------------------------------------------
  function get_angular_rate_dot(theta, dtheta) result(SMAT_DOT)

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
         & 0.0_dp, -s1*dtheta(1), c1*c2*dtheta(1) - s1*s2*dtheta(2),&
         &  0.0_dp, -c1*dtheta(1), -s1*c2*dtheta(1) - c1*s2*dtheta(2) /) )

  end function get_angular_rate_dot

  !**************************************************************
  ! create a body using the supplied state and other parameters
  !**************************************************************
  subroutine create_body(mass, re, q, qdot, alpha)

    real(dp), intent(in)        :: mass 
    type(vector), intent(in)    :: re 
    real(dp), intent(in)        :: q(NUM_STATES), qdot(NUM_STATES)
    type(body), intent(out)     :: alpha
    real(dp)                    :: theta(NUM_SPAT_DIM)

    call set_state(q, qdot, alpha)

    ! make a local copy
    theta = q(4:6)

    ! rotation and rate matrices
    alpha%C_mat     = get_rotation(theta)
    alpha%S         = get_angular_rate(theta)
    alpha%S_dot     = get_angular_rate_dot(theta, qdot(4:6))

    ! mass of the body
    alpha%mass  = mass

    ! mass moment of intertia (second moment)
    alpha%J = -mass*skew(re)*skew(re)

    ! first moment of mass
    alpha%c  = mass*re

    ! reaction force 
    alpha%fr = mass*GRAV !? check

    ! reaction torque
    alpha%gr = alpha%J*alpha%omega_dot !? check

    !   call print_body(alpha)

  end subroutine create_body


  ! ********************************************************
  ! routine that prints the state and properties of the body
  ! ********************************************************
  subroutine print_body(alpha)

    use dispmodule

    type(body):: alpha

    call disp('.....................Body.........................')
    call DISP('')
    call disp('   r          =   ', get_array(alpha%r), SEP=', ', ORIENT = 'ROW')
    call disp('   theta      =   ', get_array(alpha%theta), SEP=', ', ORIENT = 'ROW')
    call disp('   v          =   ', get_array(alpha%r), SEP=', ', ORIENT = 'ROW')
    call disp('   omega      =   ', get_array(alpha%theta), SEP=', ', ORIENT = 'ROW')
    call DISP('')
    call disp('   r_dot      =   ', get_array(alpha%r), SEP=', ', ORIENT = 'ROW')
    call disp('   theta_dot  =   ', get_array(alpha%theta), SEP=', ', ORIENT = 'ROW')
    call disp('   v_dot      =   ', get_array(alpha%r), SEP=', ', ORIENT = 'ROW')
    call disp('   omega_dot  =   ', get_array(alpha%theta), SEP=', ', ORIENT = 'ROW')
    call DISP('')
    call disp('   c          =   ', get_array(alpha%c), SEP=', ', ORIENT = 'ROW')
    call DISP('')
    call disp('   J          =   ', get_matrix(alpha%J))
    call DISP('')
    call disp('   C          =   ', get_matrix(alpha%C_mat))
    call DISP('')
    call disp('   S          =   ', get_matrix(alpha%S))
    call DISP('')
    call disp('   S_dot      =   ', get_matrix(alpha%S_dot))
    call DISP('')
    call disp('   fr         =   ', get_array(alpha%fr), SEP=', ', ORIENT = 'ROW')
    call disp('   gr         =   ', get_array(alpha%gr), SEP=', ', ORIENT = 'ROW')
    call disp('')
    call disp('..................................................')

  end subroutine print_body

end module bodies