module residual

!  use global_constants
!  use global_variables
!  use utils
  use types
  use bodies
  use joints

  implicit none

contains

  function  get_joint_residual(alpha, jnt) result(res_joint)

    type(vector) :: res_joint(NUM_JOINT_EQN)
    type(body)   :: alpha
    type(joint)   :: jnt
    

  end function get_joint_residual

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

end module residual
