!=====================================================================!
! Module that assembles the residual
!---------------------------------------------------------------------!
! The functions to call are by the name of the interface. 
!
! It is recommended to refrain from calling by name of the internal 
! functions directly as their signatures may change in future to adapt
! different requirements and scenarios.
!
!---------------------------------------------------------------------!
! This module is implemented with the following preference in mind.
!
! Prefers passing derived datatypes such as BODY and JOINT instead of 
! data in primitive form such as double arrays (Q, QDOT) because:
!
! (a) it is cleaner to pass data in and out as BEANS compared to
! primitive datatypes --that are more susceptible to programmer errors
! like passing the wrong arguments
!
! (b) future enhancements come easy when data is passed in BEAN form 
! and can be easily done without having to change the function 
! signatures. One can just add extra variables into the BODY or JOINT
! types instead of sending one or more variables as new arguments.
! For example: When extending from rigid-case to elastic case, the 
! second derivatives of the state (q_double_dot) are needed to assemble
! the residual, which may warrant a function signature change when 
! data is passed in an out in primitive form.
! 
! (c) For multi-body systems it is cleaner to keep track of bodies and 
! joints, compared to a very long Q, Q_dot, Q_double_dot arrays whose
! indices may correspond to any body or joint eqn.
!
!---------------------------------------------------------------------!
!  Example usage:
!  res(1:12) = get_residual(body)
!  res(:)    = get_residual(body_array, joint_array) (not fully func)
!---------------------------------------------------------------------!
! Inputs:
! array of bodies and joints (multi-body dynamics) OR (not fully func)
! single body (single body dynamics)
!---------------------------------------------------------------------!
! Output:
! Residual array of size NUM_STATES
!=====================================================================!
module residual

  use global_constants
  use utils
  use types
  use bodies
  use joints

  implicit none

  ! Makes all functions private by default
  private

  ! Expose only the needed functions. 
  public get_residual, get_q, get_qdot

  !*******************************************************************!
  !  Function to be called by the end-programs to get the residual
  !-------------------------------------------------------------------!
  !  Example:
  !  res(1:12) = get_residual(body)
  !  res(:)    = get_residual(body_array, joint_array) (not fully func)
  !-------------------------------------------------------------------! 
  ! This residual assembly implementation is designed assuming:
  ! (a) When a single body is supplied, the joint would not be there
  ! (b) When array of bodies is passed, it has to have some joints
  !-------------------------------------------------------------------! 
  ! The user will call this interface to get the residual array. 
  ! This interface will adapt to different types of inputs:
  !
  ! (a) It can take an array of bodies and joints for multi-body cases
  ! (b) It can take a single body too, for simple dynamics of a body     
  !-------------------------------------------------------------------! 
  ! Note: This interface does not return the NEGATIVE of the residual
  ! as one may need in the solution of linear system for update.
  !*******************************************************************!
  interface get_residual
     module procedure get_body_residual, residual_assembler
  end interface get_residual

contains
  
  !*******************************************************************!
  ! This function internally makes use of other private functions
  ! to assemble the residual and returns an ARRAY to the calling prgrm
  !-------------------------------------------------------------------!
  ! Input:  Body array and joint array
  ! Output: Full residual ARRAY
  !-------------------------------------------------------------------!
  ! Example usage:
  ! res_arr = get_residual(body_array, joint_array)
  !*******************************************************************!
  function residual_assembler(body_array, joint_array) result(res_arr)

    ! inputs
    type(body), dimension(:)             :: body_array
    type(joint), dimension(:)            :: joint_array

    ! outputs
    real(dp), dimension(TOT_NDOF)        :: res_arr

    ! get in vector form and covert to plain array form to return
    res_arr = array(residual_assembler_vec(body_array, joint_array))

  end function residual_assembler


  !*******************************************************************!
  ! This function internally makes use of other private functions
  ! to assemble the residual and returns to the calling program
  !
  ! DRAFT Implementation:
  ! Need to decide upon the order or elements (joints and bodies)
  !-------------------------------------------------------------------!
  ! Input:  Body array and joint array
  ! Output: Full residual in vector form
  !-------------------------------------------------------------------!
  ! Example usage:
  ! res_arr = array(residual_assembler_vec(body_array, joint_array))
  ! res_vec = residual_assembler_vec(body_array, joint_array)
  !*******************************************************************!
  function residual_assembler_vec(body_array, joint_array) result(res)

    ! inputs
    type(body), dimension(:)             :: body_array
    type(joint), dimension(:)            :: joint_array

    ! outputs
    type(vector), dimension(TOT_NEQN)    :: res

    ! temp local variables
    integer(sp)                          :: njnt, nbdy
    integer(sp)                          :: j, is, ie 

    !-----------------------------------------------------------------
    ! Complex joints
    !-----------------------------------------------------------------
    ! Should think how we should handle complex jonits without 
    ! affecting sparsity of the system
    !-----------------------------------------------------------------
    if (CMPLX_JOINT) then 

       return
    end if

    !-----------------------------------------------------------------
    ! Simple joints between two bodies (Need to decide order here too)
    !-----------------------------------------------------------------

    ! (a) BODY1, BODY2, JOINT (order)


    ! (b) ALL BODIES1, ALL JOINTS (will affect sparsity)

    ! first parse the bodies and extract the residual
    nbdy = size(body_array)

    ! sanity check
    if (nbdy .ne. NUM_BODIES) stop " >> Error in the number of bodies"

    is = 1
    do j = 1, nbdy
       ie = j*NUM_BODY_EQN 
       res(is:ie) = get_body_residual_vec(body_array(j))
       is = ie + 1
    end do

    ! joint residuals are put after the body residuals
    njnt = size(joint_array)

    ! sanity check
    if (njnt .ne. NUM_JOINTS) stop " >> Error in the number of joints"

    ! Note: using the same 'is' from above do loop
    do j = 1, njnt
       ie = j*NUM_JOINT_EQN 
       res(is:ie) = get_joint_residual_vec(joint_array(j))
       is = ie + 1
    end do

    stop"Incomplete Impl"

  end function residual_assembler_vec

  !*******************************************************************!
  ! Returns the residual terms coming from joint equations
  !-------------------------------------------------------------------!
  ! Input: JOINT object
  ! Output: Residual vector block in vector form
  !-------------------------------------------------------------------!
  ! Example usage:
  ! res_arr = array(get_joint_residual_vec(joint1)) OR
  ! res_vec = get_joint_residual_vec(joint1)
  !*******************************************************************!
  function  get_joint_residual_vec( jnt) result(res_jnt)

    type(joint)  :: jnt ! input joint
    type(vector) :: res_jnt(NUM_JOINT_EQN) ! residual of the joint

!    res_jnt(1) = zeroV
!    res_jnt(2) = zeroV

    stop "dummy impl"

  end function get_joint_residual_vec
  

  !*******************************************************************!
  ! Assembles the residual vector for the body from kinematics, 
  ! dynamics and elastic equations
  !-------------------------------------------------------------------!
  ! Input: BODY object
  ! Output: Residual vector block in vector form
  !-------------------------------------------------------------------!
  ! Example usage:
  ! res_arr = get_residual(body1)
  !*******************************************************************!
  function get_body_residual(alpha) result(res_dyn_arr)

    ! input
    type(body), intent(in)   :: alpha

    ! output
    real(dp), dimension(NDOF_PBODY)    :: res_dyn_arr
    real(dp), dimension(NDOF_PBODY) :: q

    res_dyn_arr = array(get_body_residual_vec(alpha))

  end function get_body_residual

  !*******************************************************************!
  ! Function that returns the state from the body
  !-------------------------------------------------------------------!
  ! Input: The body alpha
  ! Output: The state q
  !*******************************************************************!
  function get_q(alpha) result(q)
    
    type(body), intent(in)          :: alpha
    real(dp), dimension(NDOF_PBODY) :: q
    
    q=  array((/alpha%r , alpha%theta, alpha%v, alpha%omega /))

  end function get_q

  !*******************************************************************!
  ! Function that returns the time-derivative of state from the body
  !-------------------------------------------------------------------!
  ! Input: The body alpha
  ! Output: The state derivative with time qdot
  !*******************************************************************!
  function get_qdot(alpha) result(qdot)

    type(body), intent(in)          :: alpha
    real(dp), dimension(NDOF_PBODY) :: qdot

    qdot=  array((/alpha%r_dot , alpha%theta_dot, &
         &alpha%v_dot, alpha%omega_dot /))
    
  end function get_qdot

  !*******************************************************************!
  ! Actual function that assembles the residual vector for the body 
  ! from kinematics, dynamics and elastic equations.
  ! All the other functions are wrappers around this function.
  !-------------------------------------------------------------------!
  ! Input : BODY object
  ! Output: Residual vector block in vector form
  !-------------------------------------------------------------------!
  ! Example usage:
  ! res_arr = array(get_body_residual_vec(body1))
  ! res_vec = get_body_residual_vec(body1)
  !*******************************************************************!
  function get_body_residual_vec(alpha) result(res_dyn)

    ! input
    type(body), intent(in)   :: alpha

    ! output
    type(vector) :: res_dyn(NUM_BODY_EQN)

    ! create local variables (for defn of these variables look 
    ! into types.f90)
    type(vector) :: r, theta, v, omega
    type(vector) :: r_dot,  theta_dot, v_dot, omega_dot
    type(vector) :: qs, qs_dot, qs_double_dot
    type(vector) :: c, fr, gr, f 
    type(matrix) :: J, p, h, K, M, C_mat, S, S_dot
    real(dp)     :: mass

    !--------------------------------------
    ! set the values for local variables
    !--------------------------------------

    ! rigid body state variables
    r             = alpha%r
    theta         = alpha%theta
    v             = alpha%v
    omega         = alpha%omega

    ! time derivative of states
    r_dot         = alpha%r_dot
    theta_dot     = alpha%theta_dot
    v_dot         = alpha%v_dot
    omega_dot     = alpha%omega_dot

    ! elatic state variables
    qs            = alpha%qs
    qs_dot        = alpha%qs_dot
    qs_double_dot = alpha%qs_double_dot

    ! other body Attributes
    mass          = alpha%mass
    c             = alpha%c
    J             = alpha%J

    p             = alpha%p
    h             = alpha%h

    K             = alpha%K
    M             = alpha%M

    C_mat         = alpha%C_mat
    S             = alpha%S
    S_dot         = alpha%S_dot

    fr            = alpha%fr
    gr            = alpha%gr

    f             = alpha%f

    !-----------------------------------------------------------------!
    ! Now assembling the residual terms. The final vector form of 
    ! the residual can be converted into array form just by using the 
    ! array(res_vector).
    !-----------------------------------------------------------------!

    ! kineamtics eqn-1 (2 terms)
    res_dyn(1)  = C_mat*r_dot - v

    ! kinematics eqn-2 (2 terms)
    res_dyn(2)  = S*theta_dot - omega

    ! dynamics eqn-1 (7 terms)
    res_dyn(3)  = mass*v_dot - skew(c)*omega_dot +p*qs_double_dot &
         &+ skew(omega)*(mass*v - skew(c)*omega + p*qs_dot) - fr

    ! dynamics eqn 2 (8-terms)
    res_dyn(4)  = skew(c)*v_dot + J*omega_dot + h*qs_double_dot & 
         &+ skew(c)*skew(omega)*v + skew(omega)*J*omega &
         &+ skew(v)*p*qs_dot + skew(omega)*h*qs_dot -gr

    ! include elastic eqn (5 terms)
    if (ELASTIC) then
!       res_dyn(5) = trans(p)*v + trans(h)*omega_dot + M*qs_double_dot &
!            &+K*qs - f
    end if

  end function get_body_residual_vec

end module residual

!!$  !*******************************************************************
!!$  ! returns the residual terms due to the elastic equation
!!$  !*******************************************************************
!!$  function get_elastic_residual(alpha) result(res_elastic)
!!$
!!$    type(vector) :: res_elastic(NUM_ELAST_EQN)
!!$    type(body)   :: alpha
!!$
!!$    stop"dummy impl"
!!$     res_elastic(1) = trans(p)*v + trans(h)*omega_dot + M*qs_double_dot &
!!$         &+K*qs -f  ! 5 terms
!!$  end function get_elastic_residual

!!$
!!$  ! ******************************************************************
!!$  ! returns the assembled residual vector for the given body
!!$  ! Joint residual terms are implemented in the JOINT module
!!$  ! ******************************************************************
!!$  function  get_body_residual(alpha) result(res)
!!$
!!$    type(body)   :: alpha
!!$    type(vector) :: res_vector_form(NUM_GOV_EQN)
!!$    real(dp), dimension(NUM_STATES) :: res_array
!!$
!!$    res_vector_form(1:NUM_DYNAM_EQN)  = get_dynamics_residual(alpha)
!!$
!!$    if (NUM_ELAST_EQN .gt. 0 .and. elastic) then
!!$       res_body(NUM_DYNAM_EQN+1:NUM_ELAST_EQN)   =  &
!!$            &get_elastic_residual(alpha)  
!!$    end if
!!$
!!$    res_array=array(res_vector_form)
!!$
!!$  end function get_body_residual
