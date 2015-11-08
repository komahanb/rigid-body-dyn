module residual

  use global_constants
  !  use global_variables
  use utils
  use types
  use bodies
  use joints

  implicit none

  ! makes all functions private by default
  private

  ! expose only the needed functions
  public get_residual

  !---------------------------------------------------------------
  !  Function to be called by the end-programs to get the residual
  !---------------------------------------------------------------
  interface get_residual
     module procedure get_body_residual, get_joint_residual, get_residual_body_joint
  end interface get_residual

contains

  !--------------------------------------------------------------
  ! This function internally makes use of other private functions
  ! to assemble the residual vector.
  !--------------------------------------------------------------
  function get_residual_body_joint(body_array, joint_array) result(res)

    type(body), dimension(:)                 :: body_array
    type(joint), dimension(:), optional      :: joint_array
    real(dp), dimension(NUM_STATES)          :: res
    integer(sp)                              :: j, njnt, nbdy, is, ie 

    !    type(vector)                             :: get_body_residual, get_joint_residual

    ! first parse the bodies and extract the residual
    nbdy = size(body_array)
    do j = 1, nbdy
       call split(j,is,ie) ! split j index storage
       res(is:ie) = get_body_residual(body_array(j))
    end do


    ! if joints are supplied parse them and extract the residual
    ! joint residuals are put after the body residuals
    if(present(joint_array))then

       njnt = size(joint_array)
       do j = 1, njnt
          call split(nbdy+j,is,ie) ! split j index storage
          res(is:ie) = get_joint_residual(joint_array(j))
       end do

    end if

  end function get_residual_body_joint

  !---------------------------------------------------------
  ! Returns the residual terms coming from joint equations
  !---------------------------------------------------------
  function  get_joint_residual( jnt) result(res)

    type(vector) :: res_joint(NUM_JOINT_EQN)
    type(joint)  :: jnt
    real(dp), dimension(NUM_STATES)          :: res

    stop "dummy impl"

  end function get_joint_residual

  ! *********************************************************
  ! returns the assembled residual vector for the given body
  ! Joint residual terms are implemented in the JOINT module
  ! ********************************************************
  function  get_body_residual(alpha) result(res)

    type(body)   :: alpha
    type(vector) :: res_body(NUM_GOV_EQN)
    real(dp), dimension(NUM_STATES)          :: res

    res_body(1:NUM_DYNAM_EQN)  = get_dynamics_residual(alpha)

!!$    if (NUM_ELAST_EQN .gt. 0 .and. elastic) then
!!$       res_body(NUM_DYNAM_EQN+1:NUM_ELAST_EQN)   = get_elastic_residual(alpha)  
!!$    end if

    res=array(res_body)


  end function get_body_residual

  ! **********************************************************************************
  ! assembles the residual vector for the body from kinematics and dynamics equations
  ! *********************************************************************************
  function get_dynamics_residual(alpha) result(res_dyn)

    type(body), intent(in)   :: alpha
    type(vector) :: res_dyn(NUM_DYNAM_EQN)

    ! create local variables (for defn of these variables look into types.f90)
    type(vector) :: r, theta, v, omega
    type(vector) :: r_dot,  theta_dot, v_dot, omega_dot
    type(vector) :: qs, qs_dot, qs_double_dot
    type(vector) :: c, fr, gr
    type(matrix) :: J, p, h, K, M, C_mat, S, S_dot
    real(dp)     :: mass

    !--------------------------------------
    ! set the values for local variables
    !--------------------------------------

    ! rigid body state variables
    r =  alpha%r
    theta =alpha%theta
    v = alpha%v
    omega = alpha%omega

    ! time derivative of states
    r_dot =alpha%r_dot
    theta_dot = alpha%theta_dot
    v_dot = alpha%v_dot
    omega_dot = alpha%omega_dot

    ! elatic state variables
    qs = alpha%qs
    qs_dot =alpha%qs_dot
    qs_double_dot = alpha%qs_double_dot

    ! other body Attributes
    mass  =alpha%mass
    c =alpha%c
    J =alpha%J

    p  = alpha%p
    h  = alpha%h

    K  =alpha%K
    M  =alpha%M

    C_mat =alpha%C_mat
    S =alpha%S
    S_dot =alpha%S_dot

    fr = alpha%fr
    gr = alpha%gr

    !--------------------------------------
    ! now assemble the dynamics residual
    !--------------------------------------

    res_dyn(1)  = C_mat*r_dot - v

    res_dyn(2)  = S*theta_dot - omega

    res_dyn(3)  = mass*v_dot - skew(c)*omega_dot +p*qs_double_dot &
         &+ skew(omega)*(mass*v - c*omega + p*qs_dot) - fr

    res_dyn(4)  = skew(c)*v_dot + J*omega_dot + h*qs_double_dot &
         &+ skew(c)*skew(omega)*v + skew(omega)*J*omega &
         &+ skew(v)*p*qs_dot + skew(omega)*h*qs_dot &
         &+ skew(omega)*h*qs_dot -gr

    ! I might do the elastic assembly here itself

  end function get_dynamics_residual

  !*********************************************************
  ! returns the residual terms due to the elastic equation
  !******************************************************
  function get_elastic_residual(alpha) result(res_elastic)

    type(vector) :: res_elastic(NUM_ELAST_EQN)
    type(body)   :: alpha

    stop"dummy impl"

  end function get_elastic_residual

end module residual
