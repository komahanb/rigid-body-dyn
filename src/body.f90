! module that defines the configuration of a single body
module body_mod

  use types
  use global_constants
  use global_variables
  use utils

  implicit none

  ! Default make everything private
  !  private 

  ! Expose only the needed variables and functions
  !  public

contains

  ! *********************************************************
  ! returns the assembled residual vector for the given body
  ! Joint residual terms are implemented in the JOINT module
  ! ********************************************************
  function  get_body_residual(alpha) result(res_body)

    type(body)   :: alpha
    type(vector) :: res_body(NUM_GOV_EQN)

    res_body(1:NUM_DYNAM_EQN)  = get_dynamics_residual(alpha)

    if (NUM_ELAST_EQN .gt. 0) then
       res_body(NUM_DYNAM_EQN+1:NUM_ELAST_EQN)   = get_elastic_residual(alpha)  
    end if

  end function get_body_residual

  ! **********************************************************************************
  ! assembles the residual vector for the body from kinematics and dynamics equations
  ! *********************************************************************************
  function get_dynamics_residual(alpha) result(res_dyn)

    type(body)   :: alpha
    type(vector) :: res_dyn(NUM_DYNAM_EQN)

    res_dyn(1)  = alpha%C_mat*alpha%r_dot - alpha%v

    res_dyn(2)  = alpha%S*alpha%theta_dot - alpha%omega

    res_dyn(3)  = alpha%mass*alpha%v_dot - skew(alpha%c)*alpha%omega_dot +alpha%p*alpha%qs_double_dot &
         &+ skew(alpha%omega)*(alpha%mass*alpha%v - alpha%c*alpha%omega + alpha%p*alpha%qs_dot) - alpha%fr

    res_dyn(4)  = skew(alpha%c)*alpha%v_dot + alpha%J*alpha%omega_dot + alpha%h*alpha%qs_double_dot &
         &+ skew(alpha%c)*skew(alpha%omega)*alpha%v + skew(alpha%omega)*alpha%J*alpha%omega &
         &+ skew(alpha%v)*alpha%p*alpha%qs_dot + skew(alpha%omega)*alpha%h*alpha%qs_dot &
         &+ skew(alpha%omega)*alpha%h*alpha%qs_dot -alpha%gr

  end function get_dynamics_residual

  !*********************************************************
  ! returns the residual terms due to the elastic equation
  !******************************************************
  function get_elastic_residual(alpha) result(res_elastic)

    type(vector) :: res_elastic
    type(body)   :: alpha

    stop"dummy impl"

  end function get_elastic_residual

  !****************************************************
  ! Takes the body and updates/sets the state variables
  !****************************************************
  subroutine set_state(q, qdot, alpha)

    real(dp), intent(in)          :: q(NUM_STATES)
    real(dp), intent(in)          :: qdot(NUM_STATES)
    type(body),intent(inout)      :: alpha

    ! set the state
    alpha%r         = vector(q(1:3))
    alpha%theta     = vector(q(4:6))
    alpha%v         = vector(q(7:9))
    alpha%omega     = vector(q(10:12))

    ! set the time derivatives of state
    alpha%r_dot     = vector(qdot(1:3))
    alpha%theta_dot = vector(qdot(4:6))
    alpha%v_dot     = vector(qdot(7:9))
    alpha%omega_dot = vector(qdot(10:12))


    ! set the elastic state variables in the body
    !    alpha%qs        = vector(qdot(IS_QS:IE_QE))
    !    alpha%qs_dot    = vector(qdot(IS_QS_DOT:IE_QS_DOT))


  end subroutine set_state

  ! ********************************************************
  ! routine that returns the rotation matrix (euler angles)
  ! ********************************************************
  function get_rotation(theta)

    real(dp)     :: theta(NUM_SPAT_DIM)    
    type(matrix) :: get_rotation

    stop "dummy impl"

  end function get_rotation

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
  function get_angular_rate(theta)

    real(dp) :: theta(NUM_SPAT_DIM)    
    type(matrix) :: get_angular_rate

    stop "dummy impl"

  end function get_angular_rate

  ! **************************************************************************
  ! routine that returns the time derivative of ang rate matrix (euler angles)
  ! *************************************************************************
  function get_angular_rate_dot(theta)

    real(dp) :: theta(NUM_SPAT_DIM)
    type(matrix) :: get_angular_rate_dot
    
    stop "dummy impl"

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

    theta = q(4:6)

    ! rotation and rate matrices
    alpha%C_mat     = get_rotation(theta)
    alpha%S         = get_angular_rate(theta)
    alpha%S_dot     = get_angular_rate_dot(theta)

    ! mass of the body
    alpha%mass  = mass

    ! mass moment of intertia (second moment)
    alpha%J = -mass*skew(re)*skew(re)

    ! first moment of mass
    alpha%c  = mass*re

    ! reaction force 
    alpha%fr = mass*GRAV

    ! reaction torque
    alpha%gr = alpha%J*alpha%omega_dot

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

end module body_mod
